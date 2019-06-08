#!/usr/bin/env python3
## vi: tabstop=4 shiftwidth=4 softtabstop=4 expandtab
import numpy as np

from scipy import constants
from matplotlib import pyplot as plt

import adcc

from mpi4py import MPI
import veloxchem as vlx
import os
import tempfile


# Hartree to eV
eV = constants.value("Hartree energy in eV")


def plot_spectrum(energies, strengths, width=0.045):
    energies = np.array(energies).flatten()
    strengths = np.array(strengths).flatten()

    def ngauss(en, osc, x, w):
        """A normalised Gaussian"""
        fac = osc / np.sqrt(2 * np.pi * w**2)
        return fac * np.exp(-(x - en)**2 / (2 * w**2))

    xval = np.arange(np.min(np.min(energies) - 1, 0),
                     np.max(energies) + 1, 0.01)
    yval = np.zeros(xval.size)
    for en, osc in zip(energies, strengths):
        yval += ngauss(en, osc, xval, width)
    plt.plot(xval, yval)


# Run SCF in VeloxChem
with tempfile.TemporaryDirectory() as tmpdir:
    infile = os.path.join(tmpdir, "vlx.in")
    outfile = os.path.join(tmpdir, "/dev/null")

    with open(infile, "w") as fp:
        fp.write("""
                 @jobs
                 task: hf
                 @end

                 @method settings
                 basis: cc-pvdz
                 @end

                 @molecule
                 charge: 0
                 multiplicity: 1
                 units: bohr
                 xyz:
                 O 0 0 0
                 H 0 0 1.795239827225189
                 H 1.693194615993441 0 -0.599043184453037
                 @end
                 """)
    task = vlx.MpiTask([infile, outfile], MPI.COMM_WORLD)
    scfdrv = vlx.ScfRestrictedDriver(task.mpi_comm, task.ostream)
    scfdrv.conv_thresh = 1e-9
    scfdrv.compute(task.molecule, task.ao_basis, task.min_basis)
    scfdrv.task = task

print(adcc.banner())

# Run an adc2 calculation:
state = adcc.adc2(scfdrv, n_singlets=7, conv_tol=1e-8)
state = adcc.attach_state_densities(state)

#
# Get HF density matrix and nuclear dipole
#
ρ_hf_tot = np.add(*[dm for dm in scfdrv.scf_tensors['D']])

# Compute dipole integrals
dip_ao = np.array(scfdrv.scf_tensors['Mu'])
mol = scfdrv.task.molecule
nuc_charges = mol.elem_ids_to_numpy()
coords = np.vstack((mol.x_to_numpy(), mol.y_to_numpy(), mol.z_to_numpy())).T
dip_nucl = np.einsum('i,ix->x', nuc_charges, coords)

# compute nuclear dipole
# charges = mol.atom_charges()
# coords = mol.atom_coords()
# dip_nucl = np.einsum('i,ix->x', charges, coords)

#
# MP2 density correction
#
mp2dm_mo = state.ground_state.mp2_diffdm
dm_mp2_ao = mp2dm_mo.transform_to_ao_basis(state.reference_state)
ρ_mp2_tot = (dm_mp2_ao[0] + dm_mp2_ao[1]).to_ndarray() + ρ_hf_tot

#
# Compute properties
#
exc_energies = []    # Excitation energies
osc_strengths = []    # Oscillator strength

print()
print("  st  ex.ene. (au)         f     transition dipole moment (au)"
      "        state dip (au)")
for i, ampl in enumerate(state.eigenvectors):
    # Compute transition density matrix
    tdm_mo = state.ground_to_excited_tdms[i]
    tdm_ao = tdm_mo.transform_to_ao_basis(state.reference_state)
    ρ_tdm_tot = (tdm_ao[0] + tdm_ao[1]).to_ndarray()

    # Compute transition dipole moment
    tdip = np.einsum('xij,ij->x', dip_ao, ρ_tdm_tot)
    osc = 2. / 3. * np.linalg.norm(tdip)**2 * np.abs(state.eigenvalues[i])

    # Compute excited states density matrix and excited state dipole moment
    opdm_mo = state.state_diffdms[i]
    opdm_ao = opdm_mo.transform_to_ao_basis(state.reference_state)
    ρdiff_opdm_ao = (opdm_ao[0] + opdm_ao[1]).to_ndarray()
    sdip_el = np.einsum('xij,ij->x', dip_ao, ρdiff_opdm_ao + ρ_mp2_tot)
    sdip = sdip_el - dip_nucl

    # Print findings
    fmt = "{0:2d}  {1:12.8g} {2:9.3g}   [{3:9.3g}, {4:9.3g}, {5:9.3g}]"
    fmt += "   [{6:9.3g}, {7:9.3g}, {8:9.3g}]"
    print(state.kind[0], fmt.format(i, state.eigenvalues[i], osc, *tdip, *sdip))


    # Save oscillator strength and excitation energies
    osc_strengths.append(osc)
    exc_energies.append(state.eigenvalues[i])
exc_energies = np.array(exc_energies)
osc_strengths = np.array(osc_strengths)

# Plot a spectrum
plot_spectrum(exc_energies * eV, osc_strengths)
plt.xlabel("Excitation energy in eV")
plt.savefig("spectrum.pdf")