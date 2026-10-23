
#!/usr/bin/env python3

import sys

def convert(input_file, output_file, skip=1):
    if skip < 1:
        raise ValueError("skip must be >= 1")

    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        frame = 0
        written = 0

        while True:
            line = fin.readline()
            if not line:
                break

            if not line.startswith("ITEM: TIMESTEP"):
                continue

            timestep = fin.readline().strip()

            line = fin.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Expected ITEM: NUMBER OF ATOMS")
            natoms = int(fin.readline().strip())

            line = fin.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Expected ITEM: BOX BOUNDS")
            for _ in range(3):
                fin.readline()

            line = fin.readline()
            if not line.startswith("ITEM: ATOMS"):
                raise ValueError("Expected ITEM: ATOMS")

            columns = line.split()[2:]

            if all(c in columns for c in ("x", "y", "z")):
                xcol, ycol, zcol = (columns.index(c) for c in ("x", "y", "z"))
                scaled = False
            elif all(c in columns for c in ("xu", "yu", "zu")):
                xcol, ycol, zcol = (columns.index(c) for c in ("xu", "yu", "zu"))
                scaled = False
            elif all(c in columns for c in ("xs", "ys", "zs")):
                xcol, ycol, zcol = (columns.index(c) for c in ("xs", "ys", "zs"))
                scaled = True
            else:
                raise ValueError("Could not find supported coordinate columns.")

            idcol = columns.index("id") if "id" in columns else None
            typecol = columns.index("type") if "type" in columns else None
            elemcol = columns.index("element") if "element" in columns else None

            if typecol is None and elemcol is None:
                raise ValueError("Need an 'element' or 'type' column.")

            # Read every frame, but only write every skip-th frame.
            write_frame = (frame % skip == 0)
            atoms = []

            for i in range(natoms):
                fields = fin.readline().split()
                if len(fields) < len(columns):
                    raise ValueError("Unexpected end of atom data.")

                if not write_frame:
                    continue

                if elemcol is not None:
                    symbol = fields[elemcol]
                else:
                    atom_type = int(fields[typecol])
                    type_to_element = {1: "C", 2: "O"}
                    symbol = type_to_element.get(atom_type, str(atom_type))

                x = float(fields[xcol])
                y = float(fields[ycol])
                z = float(fields[zcol])
                atom_id = int(fields[idcol]) if idcol is not None else i

                atoms.append((atom_id, symbol, x, y, z))

            if write_frame:
                atoms.sort(key=lambda atom: atom[0])
                fout.write(f"{natoms}\n")
                fout.write(f"Timestep {timestep}\n")
                for _, symbol, x, y, z in atoms:
                    fout.write(f"{symbol:<3s} {x:16.8f} {y:16.8f} {z:16.8f}\n")
                written += 1

            frame += 1

    print(f"Read {frame} frames; wrote {written} frames to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        print(f"Usage: {sys.argv[0]} traj.lammpstrj output.xyz [skip]")
        sys.exit(1)

    skip = int(sys.argv[3]) if len(sys.argv) == 4 else 1
    convert(sys.argv[1], sys.argv[2], skip)
