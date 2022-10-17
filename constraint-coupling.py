import sys
import numpy as np
import MDAnalysis as mda


def compute_coupling(constraints):

    n_con = len(constraints)
    Sdiag = [np.sqrt(1/c[0].mass +1/c[1].mass) for c in constraints]
    Sdiag_inv = [1/elem for elem in Sdiag]

    # normalized constraint direction vectors
    B = np.array([c[1].position - c[0].position for c in constraints])
    B = np.array([v/np.linalg.norm(v) for v in B])


    # CONSTRAINT COUPLING MATRIX
    A = np.zeros((n_con,n_con))
    mass_factor = np.zeros((n_con,n_con))
    # Iterate through pairs of constraints
    for i,c1 in enumerate(constraints):
        for j,c2 in enumerate(constraints):

            # the common vertices (must be only 1!)
            ag = c1.atoms.intersection(c2.atoms)

            if len(ag) != 1:
                continue

            # TODO: this requirement is due to the VIS topology in https://pubs.acs.org/doi/pdf/10.1021/ct500100f
            # where a mass-less particle is bonded to a massive one, and later potentially turned into a constraint by "all-bonds"
            # The related topological element is "CONNBONDS" in the tpr.
            if ag.masses == 0.0:
                continue

            # signs
            if (c1.indices[0] == c2.indices[0]) or (c1.indices[1] == c2.indices[1]):
                sign = -1
            else:
                sign = +1

            # mass factor
            invmass = 1/ag[0].mass
            mass_factor[i,j] = sign * (invmass*Sdiag_inv[i]*Sdiag_inv[j])

            # constraint couplings
            A[i,j] = mass_factor[i,j] * np.dot(B[i],B[j])

    return A



if __name__ == '__main__':

    tpr = sys.argv[1]
    cnf = sys.argv[2]
    # constraint_file: a text file containing the index of the bond (not the particle).
    #                  This is necessary as MDA does not store bonds/constraints separately.
    #                  Also, subsets can be selected allowing for more flexible constraint-diagnostics.
    constraint_file = sys.argv[3]

    u = mda.Universe(tpr,cnf)

    constraint_list = np.loadtxt(constraint_file, dtype=np.int)
    constraints = u.bonds[constraint_list]

    # print all constraints included in the analysis
    for c in constraints:
        print (f"Using constraint between atoms {c.atoms[0].name}--{c.atoms[1].name} ({c.indices[0]}--{c.indices[1]})")

    A = compute_coupling(constraints)

    # w: eigenvalues
    w, _ = np.linalg.eig(A)

    # ESTIMATES LINCS ORDER:
    tol = 0.4**4
    eig_max = np.abs(w).max()
    power = int(np.log(tol)/np.log(eig_max))
    print ()
    print (f"Largest eigenvalue: {eig_max}")
    print (f"Estimated LINCS order: {power}")
