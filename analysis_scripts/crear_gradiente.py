#!/usr/bin/env python3

from pathlib import Path
import re
import bisect
import numpy as np
import MDAnalysis as mda


# ============================================================
# CONFIGURACIÓN
# ============================================================

# Estructura original equilibrada

path = "/media/leandro/My Passport/MCU current/alphafold-models/charmm-gui-mcu_emre/gmx_10xDelta/"
INPUT_GRO = Path(path + "nvt6.gro")

# Topología original CHARMM-GUI
INPUT_TOP = Path(path + "topol.top")

# Archivos de salida
OUTPUT_GRO = Path(path + "gradient_10x.gro")
OUTPUT_TOP = Path(path + "topol_gradient.top")

# Semilla para reproducibilidad
SEED = 20260923

# Límites de los compartimentos externos
Z_LOWER = 50.0    # Å
Z_UPPER = 260.0   # Å

# Iones a eliminar
REMOVE_CAL = 142
REMOVE_CLA = 284


# ============================================================
# FUNCIONES
# ============================================================

def parse_molecule_itp(itp_path):
    """
    Lee un .itp y devuelve:
        nombre de la molecula
        lista de nombres de atomos
    """
    lines = itp_path.read_text().splitlines()

    molecule_name = None
    atom_names = []

    section = None

    for line in lines:
        line = line.strip()

        if not line or line.startswith(";"):
            continue

        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip().lower()
            continue

        if section == "moleculetype":
            if molecule_name is None:
                fields = line.split()
                if fields:
                    molecule_name = fields[0]

        elif section == "atoms":
            fields = line.split()

            # Formato estándar:
            # nr type resnr residue atom cgnr charge mass
            if len(fields) >= 5:
                try:
                    int(fields[0])
                    atom_names.append(fields[4])
                except ValueError:
                    pass

    return molecule_name, atom_names


def parse_topology(top_path):
    """
    Lee los #include y la sección [ molecules ].
    """

    lines = top_path.read_text().splitlines()

    includes = []

    for line in lines:
        m = re.match(r'\s*#include\s+"([^"]+)"', line)

        if m:
            includes.append(m.group(1))

    atom_names_by_molecule = {}

    for inc in includes:

        inc_path = top_path.parent / inc

        if not inc_path.exists():
            continue

        mol_name, atom_names = parse_molecule_itp(inc_path)

        if mol_name is not None and atom_names:
            atom_names_by_molecule[mol_name] = atom_names

    # --------------------------------------------------------
    # [ molecules ]
    # --------------------------------------------------------

    molecule_entries = []

    in_molecules = False

    for i, line in enumerate(lines):

        stripped = line.strip()

        if stripped.startswith("[") and stripped.endswith("]"):
            section = stripped[1:-1].strip().lower()
            in_molecules = (section == "molecules")
            continue

        if not in_molecules:
            continue

        if not stripped or stripped.startswith(";"):
            continue

        fields = stripped.split()

        if len(fields) < 2:
            continue

        name = fields[0]

        try:
            count = int(fields[1])
        except ValueError:
            continue

        molecule_entries.append({
            "line_index": i,
            "name": name,
            "count": count
        })

    return lines, molecule_entries, atom_names_by_molecule


# ============================================================
# CARGAR ESTRUCTURA
# ============================================================

print("\n==============================================")
print(" CREACIÓN DEL GRADIENTE DE Ca2+")
print("==============================================\n")

print(f"Estructura original : {INPUT_GRO}")
print(f"Topología original  : {INPUT_TOP}\n")

if not INPUT_GRO.exists():
    raise FileNotFoundError(f"No existe: {INPUT_GRO}")

if not INPUT_TOP.exists():
    raise FileNotFoundError(f"No existe: {INPUT_TOP}")


u = mda.Universe(str(INPUT_GRO))

print(f"Átomos originales: {len(u.atoms)}")
print(f"Box Z: {u.dimensions[2]:.3f} Å\n")


# ============================================================
# LEER TOPOLOGÍA
# ============================================================

top_lines, molecule_entries, atom_names_by_molecule = parse_topology(
    INPUT_TOP
)

print("Molecules encontradas en topol.top:\n")

for entry in molecule_entries:
    name = entry["name"]
    count = entry["count"]

    if name in atom_names_by_molecule:
        natoms = len(atom_names_by_molecule[name])
    else:
        natoms = "?"

    print(f"{name:8s} {count:8d}   atoms/molecule = {natoms}")


# ============================================================
# CONSTRUIR LOS BLOQUES MOLECULARES
# ============================================================

print("\nConstruyendo bloques moleculares...")

blocks = []
atom_cursor = 0

for entry in molecule_entries:

    name = entry["name"]
    count = entry["count"]

    if name not in atom_names_by_molecule:
        raise RuntimeError(
            f"No se pudo determinar el número de átomos de {name}"
        )

    atom_names = atom_names_by_molecule[name]
    atoms_per_molecule = len(atom_names)

    total_atoms = count * atoms_per_molecule

    start = atom_cursor
    end = atom_cursor + total_atoms

    blocks.append({
        "name": name,
        "count": count,
        "atoms_per_molecule": atoms_per_molecule,
        "start": start,
        "end": end,
        "line_index": entry["line_index"]
    })

    atom_cursor = end


if atom_cursor != len(u.atoms):

    raise RuntimeError(
        "\nERROR:\n"
        f"La topología predice {atom_cursor} átomos,\n"
        f"pero el GRO contiene {len(u.atoms)}.\n\n"
        "Esto indica que el orden del GRO y topol.top "
        "ya no corresponden.\n"
    )


print(f"\n✓ Topología y GRO contienen {atom_cursor} átomos.")
print("✓ El orden molecular original es consistente.")


# ============================================================
# COMPROBAR NOMBRES DE ÁTOMOS
# ============================================================

print("\nComprobando nombres de átomos...")

expected_atom_names = []

for block in blocks:

    names = atom_names_by_molecule[block["name"]]

    for _ in range(block["count"]):
        expected_atom_names.extend(names)


actual_atom_names = list(u.atoms.names)

if len(expected_atom_names) != len(actual_atom_names):
    raise RuntimeError("Error en el número de nombres de átomos.")

mismatch = []

for i, (expected, actual) in enumerate(
        zip(expected_atom_names, actual_atom_names)):

    if expected != actual:
        mismatch.append((i + 1, expected, actual))

        if len(mismatch) >= 10:
            break


if mismatch:

    print("\nPrimeros mismatches:")

    for atom, expected, actual in mismatch:
        print(
            f"Átomo {atom}: "
            f"topología={expected} "
            f"GRO={actual}"
        )

    raise RuntimeError(
        "\nEl GRO original ya no coincide con topol.top.\n"
        "No se modificará ningún archivo."
    )

print("✓ Los nombres de átomos coinciden con la topología.")


# ============================================================
# SELECCIÓN DE IONES EXTERNOS
# ============================================================

z = u.atoms.positions[:, 2]

external_mask = (
    (z < Z_LOWER) |
    (z > Z_UPPER)
)


cal = u.select_atoms("resname CAL")
cla = u.select_atoms("resname CLA")

cal_external = cal[external_mask[cal.indices]]
cla_external = cla[external_mask[cla.indices]]

print("\nIones externos disponibles:")
print(f"CAL externo: {len(cal_external)}")
print(f"CLA externo: {len(cla_external)}")


if len(cal_external) < REMOVE_CAL:
    raise RuntimeError("No hay suficientes CAL externos.")

if len(cla_external) < REMOVE_CLA:
    raise RuntimeError("No hay suficientes CLA externos.")


# ============================================================
# SELECCIÓN ALEATORIA REPRODUCIBLE
# ============================================================

rng = np.random.default_rng(SEED)

remove_cal = rng.choice(
    cal_external.indices,
    size=REMOVE_CAL,
    replace=False
)

remove_cla = rng.choice(
    cla_external.indices,
    size=REMOVE_CLA,
    replace=False
)

remove_indices = np.concatenate([
    remove_cal,
    remove_cla
])

remove_indices = np.sort(remove_indices)


print("\nEliminación:")
print(f"CAL: {len(remove_cal)}")
print(f"CLA: {len(remove_cla)}")
print(f"Total: {len(remove_indices)}")


# ============================================================
# VERIFICACIÓN DE CARGA
# ============================================================

removed_charge = (
    len(remove_cal) * 2
    - len(remove_cla)
)

print(f"\nCarga eliminada: {removed_charge:+d} e")

if removed_charge != 0:
    raise RuntimeError(
        "La eliminación de iones no conserva la neutralidad."
    )

print("✓ La eliminación conserva exactamente la carga.")


# ============================================================
# DETERMINAR QUÉ BLOQUES CAMBIAN
# ============================================================

block_ends = [b["end"] for b in blocks]


removed_by_block = {}

for idx in remove_indices:

    block_idx = bisect.bisect_right(block_ends, idx)

    if block_idx >= len(blocks):
        raise RuntimeError(
            f"No se pudo localizar el átomo {idx}"
        )

    block = blocks[block_idx]

    removed_by_block.setdefault(
        block_idx,
        {"CAL": 0, "CLA": 0}
    )

    atom_name = u.atoms[idx].resname

    if atom_name not in ("CAL", "CLA"):
        raise RuntimeError(
            f"Se intentó eliminar un residuo inesperado: {atom_name}"
        )

    removed_by_block[block_idx][atom_name] += 1


# ============================================================
# NUEVOS COUNTS DE TOPOLOGÍA
# ============================================================

new_counts = []

print("\nCambios en [ molecules ]:\n")

for i, block in enumerate(blocks):

    old_count = block["count"]

    removed_cal = removed_by_block.get(i, {}).get("CAL", 0)
    removed_cla = removed_by_block.get(i, {}).get("CLA", 0)

    removed = removed_cal + removed_cla

    # CAL y CLA son moléculas de 1 átomo
    new_count = old_count - removed

    if new_count < 0:
        raise RuntimeError(
            f"El número de {block['name']} sería negativo."
        )

    new_counts.append(new_count)

    if new_count != old_count:

        print(
            f"{block['name']:8s}: "
            f"{old_count:6d} -> {new_count:6d}"
        )


# ============================================================
# CREAR NUEVO GRO
# ============================================================

keep = np.ones(len(u.atoms), dtype=bool)

keep[remove_indices] = False

gradient = u.atoms[keep]

print(
    f"\nÁtomos después de eliminar iones: "
    f"{len(gradient)}"
)

gradient.write(str(OUTPUT_GRO))

print(f"✓ Escrito: {OUTPUT_GRO}")


# ============================================================
# CREAR NUEVA TOPOLOGÍA
# ============================================================

new_top_lines = list(top_lines)

for block, new_count in zip(blocks, new_counts):

    line_idx = block["line_index"]

    original_line = new_top_lines[line_idx]

    # Conserva comentarios si los hubiera
    match = re.match(
        r'^(\s*' +
        re.escape(block["name"]) +
        r'\s+)\d+(\s*(?:;.*)?)$',
        original_line.rstrip("\n")
    )

    if not match:
        raise RuntimeError(
            f"No se pudo modificar la línea de {block['name']} "
            f"en topol.top."
        )

    new_top_lines[line_idx] = (
        f"{match.group(1)}"
        f"{new_count}"
        f"{match.group(2)}\n"
    )


OUTPUT_TOP.write_text(
    "\n".join(new_top_lines) + "\n"
)

print(f"✓ Escrito: {OUTPUT_TOP}")


# ============================================================
# VALIDACIÓN FINAL
# ============================================================

expected_final_atoms = 0

for block, new_count in zip(blocks, new_counts):

    expected_final_atoms += (
        new_count *
        block["atoms_per_molecule"]
    )


print("\n==============================================")
print(" VALIDACIÓN FINAL")
print("==============================================")

print(
    f"Átomos esperados por topología : "
    f"{expected_final_atoms}"
)

print(
    f"Átomos presentes en GRO        : "
    f"{len(gradient)}"
)


if expected_final_atoms != len(gradient):

    raise RuntimeError(
        "\nERROR FINAL: GRO y topología no tienen "
        "el mismo número de átomos."
    )


print("✓ Número de átomos consistente.")


# Conteo final de iones
final_cal = np.sum(
    gradient.resnames == "CAL"
)

final_cla = np.sum(
    gradient.resnames == "CLA"
)

final_pot = np.sum(
    gradient.resnames == "POT"
)

print("\nIones finales:")
print(f"CAL : {final_cal}")
print(f"CLA : {final_cla}")
print(f"POT : {final_pot}")

print("\n==============================================")
print(" GRADIENTE CREADO CORRECTAMENTE")
print("==============================================")
print(f"\nGRO      : {OUTPUT_GRO}")
print(f"Topology : {OUTPUT_TOP}")
print()