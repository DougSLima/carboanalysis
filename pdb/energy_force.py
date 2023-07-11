# -*- coding: utf-8 -*-
"""
Computing Energy and Force Using Models Inside Model Zoo
========================================================

TorchANI has a model zoo trained by NeuroChem. These models are shipped with
TorchANI and can be used directly.
"""

###############################################################################
# To begin with, let's first import the modules we will use:
import torch
import torchani
import pandas as pd
import numpy as np


###############################################################################
# Let's now manually specify the device we want TorchANI to run:
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

###############################################################################
# Let's now load the built-in ANI-1ccx models. The builtin ANI-1ccx contains 8
# models trained with diffrent initialization. Predicting the energy and force
# using the average of the 8 models outperform using a single model, so it is
# always recommended to use an ensemble, unless the speed of computation is an
# issue in your application.
#
# The ``periodic_table_index`` arguments tells TorchANI to use element index
# in periodic table to index species. If not specified, you need to use
# 0, 1, 2, 3, ... to index species
model = torchani.models.ANI2x(periodic_table_index=True).to(device)

###############################################################################
# Now let's define the coordinate and species. If you just want to compute the
# energy and force for a single structure like in this example, you need to
# make the coordinate tensor has shape ``(1, Na, 3)`` and species has shape
# ``(1, Na)``, where ``Na`` is the number of atoms in the molecule, the
# preceding ``1`` in the shape is here to support batch processing like in
# training. If you have ``N`` different structures to compute, then make it
# ``N``.
#
# .. note:: The coordinates are in Angstrom, and the energies you get are in Hartree
coordinates_1 = torch.tensor([[[11.890, 13.200, 13.900],
                             [12.070, 12.640, 15.180],
                             [12.460, 12.280, 12.830],
                             [11.880, 10.970, 12.870],
                             [13.980, 12.220, 12.910],
                             [14.520, 11.400, 11.880],
                             [14.550, 13.620, 12.850],
                             [15.960, 13.600, 12.970],
                             [13.920, 14.480, 13.940],
                             [12.490, 14.500, 13.790],
                             [14.410, 15.910, 13.920],
                             [13.810, 16.690, 14.950],
                             [11.630, 13.220, 15.800],
                             [11.870, 10.670, 13.780],
                             [13.930, 10.650, 11.770],
                             [16.280, 12.880, 12.410],
                             [12.860, 16.700, 14.810],
                             [10.800, 13.310, 13.720],
                             [12.200, 12.720, 11.840],
                             [14.260, 11.760, 13.880],
                             [14.290, 14.070, 11.870],
                             [14.200, 14.050, 14.930],
                             [15.510, 15.940, 14.060],
                             [14.190, 16.370, 12.930]]],
                           requires_grad=True, device=device)

# coordinates_1 = coordinates_1.detach().numpy()

# In periodic table, C = 6 and H = 1

species_1 = torch.tensor([[6, 8, 6, 8, 6, 8, 6, 8, 6, 8, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], device=device)



###############################################################################
# Now let's compute energy and force:
energy_1_hartree = model((species_1, coordinates_1)).energies
energy_1_kjmol = torchani.units.hartree2kjoulemol(energy_1_hartree)

derivative_1 = torch.autograd.grad(energy_1_hartree.sum(), coordinates_1)[0]
force_1 = -derivative_1

###############################################################################
# And print to see the result:
print('Energy 1.pdb (Hartree):', energy_1_hartree.item())
print('Energy 1.pdb (kJ/mol): ', energy_1_kjmol.item())
#print('Force 1.pdb:', force_1.squeeze())

###############################################################################
# you can also get the atomic energies (WARNING: these have no physical
# meaning) by calling:
_, atomic_energies_1 = model.atomic_energies((species_1, coordinates_1))

###############################################################################
# this gives you the average (shifted) energies over all models of the ensemble by default,
# with the same shape as the coordinates. Dummy atoms, if present, will have an
# energy of zero
#print('Average Atomic energies, for species 6 1 1 1 1', atomic_energies_1)

###############################################################################
# you can also access model specific atomic energies
_, atomic_energies = model.atomic_energies((species_1, coordinates_1), average=False)
#print('Atomic energies of first model, for species 6 1 1 1 1', atomic_energies[0, :, :])

###############################################################################
# Same thing for 2.pdb
coordinates_2 = torch.tensor([[[12.140, 13.250, 14.090],
                             [11.080, 12.740, 14.780],
                             [12.700, 12.220, 13.070],
                             [13.080, 10.930, 13.610],
                             [13.840, 12.800, 12.220],
                             [13.290, 13.820, 11.380],
                             [14.940, 13.370, 13.150],
                             [15.730, 12.320, 13.730],
                             [14.430, 14.290, 14.310],
                             [13.230, 13.730, 14.920],
                             [14.250, 15.810, 13.990],
                             [13.390, 16.130, 12.910],
                             [11.680, 12.110, 15.210],
                             [12.820, 10.260, 12.980],
                             [13.470, 14.660, 11.800],
                             [15.170, 11.710, 14.210],
                             [12.700, 16.740, 13.170],
                             [11.760, 14.120, 13.510],
                             [11.890, 12.000, 12.330],
                             [14.210, 12.000, 11.520],
                             [15.680, 13.940, 12.530],
                             [15.190, 14.100, 15.090],
                             [13.990, 16.360, 14.910],
                             [15.320, 16.140, 13.890]]],
                           requires_grad=True, device=device)

# coordinates_2 = coordinates_2.detach().numpy()

# In periodic table, C = 6 and H = 1
species_2 = torch.tensor([[6, 8, 6, 8, 6, 8, 6, 8, 6, 8, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], device=device)

###############################################################################
# Now let's compute energy and force:
energy_2_hartree = model((species_2, coordinates_2)).energies
energy_2_kjmol = torchani.units.hartree2kjoulemol(energy_2_hartree)

print('Energy 2.pdb (Hartree):', energy_2_hartree.item())
print('Energy 2.pdb (kJ/mol): ', energy_2_kjmol.item())

derivative_2 = torch.autograd.grad(energy_2_hartree.sum(), coordinates_2)[0]
force_2 = -derivative_2

###############################################################################
# And print to see the result:
#print('Energy 2.pdb:', energy_2.item())
#print('Force 2.pdb:', force_2.squeeze())

###############################################################################
# you can also get the atomic energies (WARNING: these have no physical
# meaning) by calling:
_, atomic_energies_2 = model.atomic_energies((species_2, coordinates_2))

###############################################################################
# this gives you the average (shifted) energies over all models of the ensemble by default,
# with the same shape as the coordinates. Dummy atoms, if present, will have an
# energy of zero
#print('Average Atomic energies, for species 6 1 1 1 1', atomic_energies_2)

###############################################################################
# you can also access model specific atomic energies
_, atomic_energies = model.atomic_energies((species_2, coordinates_2), average=False)
#print('Atomic energies of first model, for species 6 1 1 1 1', atomic_energies[0, :, :])

###############################################################################
# Same thing for 3.pdb
coordinates_3 = torch.tensor([[[11.810, 12.920, 13.290],
                             [10.700, 12.720, 14.090],
                             [12.790, 11.810, 13.360],
                             [13.170, 11.440, 14.700],
                             [14.070, 12.230, 12.670],
                             [13.970, 12.030, 11.270],
                             [14.430, 13.730, 12.940],
                             [15.840, 13.880, 13.060],
                             [13.600, 14.340, 14.040],
                             [12.270, 14.260, 13.540],
                             [13.860, 15.780, 14.310],
                             [15.260, 16.030, 14.720],
                             [10.380, 13.420, 14.650],
                             [12.360, 11.090, 15.090],
                             [14.600, 12.620, 10.840],
                             [16.040, 14.260, 13.910],
                             [15.430, 16.960, 14.890],
                             [11.360, 13.030, 12.280],
                             [12.490, 10.880, 12.850],
                             [14.860, 11.550, 13.040],
                             [14.140, 14.270, 12.010],
                             [13.750, 13.760, 14.960],
                             [13.700, 16.270, 13.330],
                             [13.090, 16.080, 15.050]]],
                           requires_grad=True, device=device)


# In periodic table, C = 6 and H = 1
species_3 = torch.tensor([[6, 8, 6, 8, 6, 8, 6, 8, 6, 8, 6, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]], device=device)

###############################################################################
# Now let's compute energy and force:
energy_3_hartree = model((species_3, coordinates_3)).energies
energy_3_kjmol = torchani.units.hartree2kjoulemol(energy_3_hartree)

print('Energy 3.pdb (Hartree):', energy_3_hartree.item())
print('Energy 3.pdb (kJ/mol): ', energy_3_kjmol.item())

derivative_3 = torch.autograd.grad(energy_3_hartree.sum(), coordinates_3)[0]
force_3 = -derivative_3

###############################################################################
# And print to see the result:
#print('Energy 3.pdb:', energy_3.item())
#print('Force 3.pdb:', force_3.squeeze())

###############################################################################
# you can also get the atomic energies (WARNING: these have no physical
# meaning) by calling:
_, atomic_energies_3 = model.atomic_energies((species_3, coordinates_3))

###############################################################################
# this gives you the average (shifted) energies over all models of the ensemble by default,
# with the same shape as the coordinates. Dummy atoms, if present, will have an
# energy of zero
#print('Average Atomic energies, for species 6 1 1 1 1', atomic_energies_3)

###############################################################################
# you can also access model specific atomic energies
_, atomic_energies = model.atomic_energies((species_3, coordinates_3), average=False)
#print('Atomic energies of first model, for species 6 1 1 1 1', atomic_energies[0, :, :])

entry_names = ['1.pdb', '2.pdb', '3.pdb']
energies_hartree = [energy_1_hartree.item(), energy_2_hartree.item(), energy_3_hartree.item()]
energies_kjmol = [energy_1_kjmol.item(), energy_2_kjmol.item(), energy_3_kjmol.item()]

dict = {"entry_name": entry_names, "energies_hartree": energies_hartree, "energies_kjmol": energies_kjmol}
df = pd.DataFrame(data = dict)

df.to_csv(path_or_buf="/home/douglas/carboanalysis/carboanalysis/pdb/dataframes/torchani/energies.csv", mode='w', index=False, header=['entry_name', 'energy_hartree', 'energy_kJ/mol'], sep=";")