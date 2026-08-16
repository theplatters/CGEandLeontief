"""Clean parse of Results.txt — proper comparison."""
import re

lines = open('/workspace/BFrep/(3)BeyondHulten/Replication Files/GDP Simulatin -- 88 Sector/Robustness/Results.txt').readlines()

current_mode = 'Full Reallocation'
current_step = 1

def find_moments(mode, step, sig, eps, th):
    """Find the row matching given parameters."""
    curr_mode = ''
    curr_step = ''
    for line in lines:
        line = line.strip()
        if 'Reallocation' in line:
            curr_mode = line
            continue
        if line.startswith('Step size'):
            curr_step = line.split()[-1]
            continue
        if curr_mode != mode or curr_step != str(step):
            continue
        m = re.match(r'\(([\d.]+),\s*([\d.]+),\s*([\d.]+)\)\s*\*\s*(-?[\d.]+)\s*\*\s*([\d.]+)\s*\*\s*(-?[\d.]+)\s*\*\s*(-?[\d.]+)\s*\*\s*([\d.eE+-]+)', line)
        if m:
            s, e, t = float(m.group(1)), float(m.group(2)), float(m.group(3))
            if abs(s-sig)<0.01 and abs(e-eps)<0.01 and abs(t-th)<0.001:
                return (float(m.group(4)), float(m.group(5)), float(m.group(6)), float(m.group(7)), float(m.group(8)))
    return None

print("Paper benchmark (σ=0.9, ε=0.5, θ=0.001), No Reallocation, step=1 (annual):")
p = find_moments('No Reallocation', 1, 0.9, 0.5, 0.001)
j = (-0.003514, 0.011225, -0.1549, 0.1189, None)  # our Julia values
if p:
    print(f"  MATLAB: mean={p[0]*100:.3f}%  std={p[1]*100:.3f}%  skew={p[2]:.4f}  exkurt={p[3]:.4f}  ampl={p[4]:.4f}")
    print(f"  Julia:  mean={j[0]*100:.3f}%  std={j[1]*100:.3f}%  skew={j[2]:.4f}  exkurt={j[3]:.4f}")
    print(f"  DIFF:   mean={abs(p[0]-j[0])*100:.4f}pp  std={abs(p[1]-j[1])*100:.4f}pp  skew={abs(p[2]-j[2]):.4f}  exkurt={abs(p[3]-j[3]):.4f}")

print()
print("Paper benchmark (σ=0.9, ε=0.5, θ=0.001), No Reallocation, step=4:")
p4 = find_moments('No Reallocation', 4, 0.9, 0.5, 0.001)
if p4:
    print(f"  MATLAB: mean={p4[0]*100:.3f}%  std={p4[1]*100:.3f}%  skew={p4[2]:.4f}  exkurt={p4[3]:.4f}  ampl={p4[4]:.4f}")

print()
print("Paper benchmark (σ=0.9, ε=0.5, θ=0.001), Full Reallocation, step=4:")
pr4 = find_moments('Full Reallocation', 4, 0.9, 0.5, 0.001)
if pr4:
    print(f"  MATLAB: mean={pr4[0]*100:.3f}%  std={pr4[1]*100:.3f}%  skew={pr4[2]:.4f}  exkurt={pr4[3]:.4f}  ampl={pr4[4]:.4f}")
    # Our realloc MC used Cov4 (4-year)
    jr = (-0.009718, 0.024797, -0.1972, 0.0336, None)
    print(f"  Julia:  mean={jr[0]*100:.3f}%  std={jr[1]*100:.3f}%  skew={jr[2]:.4f}  exkurt={jr[3]:.4f}")
    print(f"  DIFF:   mean={abs(pr4[0]-jr[0])*100:.4f}pp  std={abs(pr4[1]-jr[1])*100:.4f}pp")

print()
print("Paper benchmark (σ=0.9, ε=0.5, θ=0.001), Full Reallocation, step=1:")
pr1 = find_moments('Full Reallocation', 1, 0.9, 0.5, 0.001)
if pr1:
    print(f"  MATLAB: mean={pr1[0]*100:.3f}%  std={pr1[1]*100:.3f}%  skew={pr1[2]:.4f}  exkurt={pr1[3]:.4f}  ampl={pr1[4]:.4f}")