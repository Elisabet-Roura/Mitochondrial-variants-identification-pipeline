#!/usr/bin/env python3
import os, subprocess, sys, logging, shutil

coverages   = [1, 100, 5000]
proportions = (0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)

out_root = "/media/gencardio/KINGSTON/Mitochondrial_simulation_sets_reference/"

def run(cmd):
    print(">>", cmd)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        logging.error(p.stderr.strip())
        sys.exit(f"[ERROR] fallo ejecutando: {cmd}")
    return p

def compressor():
    """Devuelve 'pigz -c' si existe pigz; si no, 'gzip -c'."""
    return "pigz -c" if shutil.which("pigz") else "gzip -c"

def sample_pair_by_prop(R1, R2, prop, seed, out_prefix):

    sub_R1  = out_prefix + ".R1.fastq"
    sub_R2  = out_prefix + ".R2.fastq"
    ids_r1  = out_prefix + ".ids.txt"
    ids_r2  = out_prefix + ".ids_r2.txt"

   
    if not os.path.exists(sub_R1):
        run(f"seqkit sample -p {prop} -s {seed} {R1} > {sub_R1}")

    if not os.path.exists(ids_r1):
        run(f"seqkit seq -n {sub_R1} > {ids_r1}")

    if not os.path.exists(ids_r2):
        with open(ids_r1, "r") as fin, open(ids_r2, "w") as fout:
            for line in fin:
                x = line.rstrip("\n")
                if x.endswith("/1"):
                    x = x[:-2] + "/2"
                fout.write(x + "\n")

    if not os.path.exists(sub_R2):
        run(f"seqkit grep -n -f {ids_r2} {R2} > {sub_R2}")

    return sub_R1, sub_R2

def main():
    os.makedirs(out_root, exist_ok=True)
    comp = compressor()  # 'pigz -c' o 'gzip -c'

    for cov in coverages:
        # Rutas de entrada
        s1_R1 = f"/media/gencardio/KINGSTON/Mitochondrial_simulation/RB{cov}_R1.fastq.gz"
        s1_R2 = f"/media/gencardio/KINGSTON/Mitochondrial_simulation/RB{cov}_R2.fastq.gz"
        s2_R1 = f"/media/gencardio/KINGSTON/Mitochondrial_simulation_ref/RB{cov}_R1.fastq.gz"
        s2_R2 = f"/media/gencardio/KINGSTON/Mitochondrial_simulation_ref/RB{cov}_R2.fastq.gz"

        for path in (s1_R1, s1_R2, s2_R1, s2_R2):
            if not os.path.exists(path):
                sys.exit(f"[ERROR] No existe: {path}")

        cov_dir   = os.path.join(out_root, f"RB{cov}")
        sets_dir  = os.path.join(cov_dir, "sets")
        mixed_dir = os.path.join(cov_dir, "mixed")
        final_dir = os.path.join(cov_dir, "final")
        for d in (sets_dir, mixed_dir, final_dir):
            os.makedirs(d, exist_ok=True)

        for p in proportions:
            q   = 1 - p
            tag = f"{p:.2f}".replace("0.", "")
            tag2 = f"{q:.2f}".replace("0.", "")

            # Subsample
            A_prefix = os.path.join(sets_dir, f"setA-p{tag}")
            B_prefix = os.path.join(sets_dir, f"setB-p{tag2}")
            A_R1, A_R2 = sample_pair_by_prop(s1_R1, s1_R2, p, 42, A_prefix)
            B_R1, B_R2 = sample_pair_by_prop(s2_R1, s2_R2, q, 42, B_prefix)

            # Mix
            mixed_R1 = os.path.join(mixed_dir, f"mixed-p{tag}-R1.fastq")
            mixed_R2 = os.path.join(mixed_dir, f"mixed-p{tag}-R2.fastq")
            if not os.path.exists(mixed_R1):
                run(f"cat {A_R1} {B_R1} > {mixed_R1}")
            if not os.path.exists(mixed_R2):
                run(f"cat {A_R2} {B_R2} > {mixed_R2}")

            # Compress to .fastq.gz
            final_R1 = os.path.join(final_dir, f"RB{cov}p{tag}_R1.fastq.gz")
            final_R2 = os.path.join(final_dir, f"RB{cov}p{tag}_R2.fastq.gz")

            if not os.path.exists(final_R1):
                run(f"{comp} {mixed_R1} > {final_R1}")
            if not os.path.exists(final_R2):
                run(f"{comp} {mixed_R2} > {final_R2}")

            print(f"[OK] {cov}x p={p}: {final_R1}, {final_R2}")

if __name__ == "__main__":
    main()
