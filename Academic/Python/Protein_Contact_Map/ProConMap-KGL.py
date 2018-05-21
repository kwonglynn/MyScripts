#!/usr/bin/env python
# encoding: utf-8
"""

Author:
    Kian Ho <hui.kian.ho@gmail.com>

Usage:
    pconpy.py cmap <dist> -p <pdb> -o <file> [options]
    pconpy.py dmap -p <pdb> -o <file> [options]
    pconpy.py hbmap -p <pdb> -o <file> [options]

Options:
    -p, --pdb <pdb>             The PDB file.
    -c, --chains <chain-ids>    Comma-separated list of chain identifiers
                                (defaults to the first chain).
    -o, --output <file>         Save the plot to a file. The file format is
                                determined by the file extension.
    -m, --measure <measure>     The inter-residue distance measure [default: CA].
    -M, --mask-thresh <dist>    Hide the distances below a given threshold (in
                                angstroms).
    --plaintext                 Generate a plaintext distance/contact matrix
                                and write to stdout (recommended for
                                piping into other CLI programs).
    --asymmetric                Display the plot only in the upper-triangle.

    --title TITLE               The title of the plot (optional).
    --xlabel <label>            X-axis label [default: Residue index].
    --ylabel <label>            Y-axis label [default: Residue index].

    --font-family <font>        Font family (via matplotlib) [default: sans].
    --font-size <size>          Font size in points [default: 10].

    --width-inches <width>      Width of the plot in inches [default: 6.0].
    --height-inches <height>    Height of the plot in inches [default: 6.0].
    --dpi <dpi>                 Set the plot DPI [default: 80]

    --greyscale                 Generate a greyscale distance map.
    --no-colorbar               Hide the color bar on distance maps.
    --transparent               Set the background to transparent.
    --show-frame

    -v, --verbose               Verbose mode.

Distance measures:
    "CA" -- Conventional CA-CA distance, this is the default distance measure.
    "CB" -- The CB-CB distance.
    "cmass" -- The distance between the residue centers of mass.
    "sccmass" -- The distance between the sidechain centers of mass
    "minvdw" -- The minimum distance between the VDW radii of each residue.

"""

import matplotlib as mpl

mpl.use('Agg') # Agg is a non-interactive backend, 
               # meaning it won't display on the screen, only save to files.

import os
import sys
import numpy
import pylab

from itertools import combinations_with_replacement
from docopt import docopt

import Bio.PDB

# The atom names of the backbone and sidechain atoms are based on those defined
# in the `ProDy` module.
BACKBONE_ATOMS = set(["CA", "C", "O", "N"])
BACKBONE_FULL_ATOMS = set(['CA', 'C', 'O', 'N', 'H', 'H1', 'H2', 'H3', 'OXT'])

# VDW radii values
#
# NOTE: 
#   - The first character in the PDB atom name is the element of the atom.
#
VDW_RADII = { "N" : 1.55, "C" : 1.70,"H" : 1.20, "O" : 1.52, "S": 1.85 }

#
# BioPython enhancements
#

def is_backbone(atom):
    return atom.get_id() in BACKBONE_ATOMS

def is_sidechain(atom):
    return atom.get_id() not in BACKBONE_FULL_ATOMS

def get_backbone_atoms(res):
    return filter(lambda atom : is_backbone(atom), res.get_iterator())

def get_sidechain_atoms(res, infer_CB=False):
    """Get the sidechain atoms of a residue.

    Args:
        res: A Bio.PDB.Residue object.
        infer_CB: Set to True to infer the CB atom if it didn't already exist
            (optional).

    Returns:
        A list of Bio.PDB.Residue objects.

    """
    atoms = filter(lambda atom : is_sidechain(atom), res.get_iterator())

    if (not atoms):
        if infer_CB:
            CB_coord = get_atom_coord(res, "CB")
            CB_atom = Bio.PDB.Atom.Atom("CB", CB_coord,
                        None, None, None, None, None, "C")
            atoms = [CB_atom]
        else:
            # Use the CA atom if there are no sidechain atoms.
            assert("CA" in res)
            atoms = [res["CA"]]

    return atoms


def get_residues(pdb_fn, chain_ids=None, model_num=0):
    """Build a simple list of residues from a single chain of a PDB file.

    Args:
        pdb_fn: The path to a PDB file.
        chain_ids: A list of single-character chain identifiers.
        model_num: The model number in the PDB file to use (optional)

    Returns:
        A list of Bio.PDB.Residue objects.

    """

    pdb_id = os.path.splitext(os.path.basename(pdb_fn))[0]

    parser = Bio.PDB.PDBParser(pdb_id, pdb_fn)
    struct = parser.get_structure(pdb_id, pdb_fn)
    model = struct[model_num]

    if chain_ids is None:
        # get residues from every chain.
        chains = model.get_list()
    else:
        chains = [ model[ch_id] for ch_id in chain_ids ]

    residues = []

    for ch in chains:
        for res in filter(lambda r : Bio.PDB.is_aa(r), ch.get_residues()):
            if not Bio.PDB.is_aa(res, standard=True):
                sys.stderr.write("WARNING: non-standard AA at %r%s"
                        % (res.get_id(), os.linesep))
            residues.append(res)

    return residues

def get_atom_coord(res, atom_name, verbose=False):
    """Get the atomic coordinate of a single atom in a residue. This function wraps the
    ``Bio.PDB.Residue.get_coord()`` function to infer CB coordinates if required.

    Args:
        res: A Bio.PDB.Residue object.
        atom_name: The name of the atom (e.g. "CA" for alpha-carbon).
        verbose: Display diagnostic messages to stderr (optional).

    Returns:
        The coordinate of the specified atom.

    """

    try:
        coord = res[atom_name].get_coord()
    except KeyError:
        if atom_name != "CB":
            # Return the first/only available atom.
            atom = res.child_dict.values()[0]
            sys.stderr.write(
                "WARNING: {} atom not found in {}".format(atom_name, res)
                    + os.linesep
            )
            return atom.get_coord()

        if verbose:
            sys.stderr.write(
                "WARNING: computing virtual {} coordinate.".format(atom_name)
                    + os.linesep)

        assert("N" in res)
        assert("CA" in res)
        assert("C" in res)

        # Infer the CB atom position described in http://goo.gl/OaNjxe
        #
        # NOTE:
        # These are Bio.PDB.Vector objects and _not_ numpy arrays.
        N = res["N"].get_vector()
        CA = res["CA"].get_vector()
        C = res["C"].get_vector()

        CA_N = N - CA
        CA_C = C - CA

        rot_mat = Bio.PDB.rotaxis(-numpy.pi * 120.0 / 180.0, CA_C)
        coord = (CA + CA_N.left_multiply(rot_mat)).get_array()

    return coord

def is_parallel(r1, r2):
    """Determine if the N-to-C directions of two residues are parallel to each other.

    """
    v1 = get_atom_coord(r1, "C") - get_atom_coord(r1, "N")
    v2 = get_atom_coord(r2, "C") - get_atom_coord(r2, "N")
    return v1.dot(v2) > 0
#
# Plotting
#

def px2pt(p):
    """Convert pixels to points.

    """
    return p * 72. / 96.


def init_spines(hidden=[]):
    """Initialise the plot frame, hiding the selected spines.

    Args:
        hidden: A list of spine names to hide. For example, set hidden
            to ["top", "right"] to hide both the top and right axes borders from
            the plot. All spines will be hidden if hidden is an empty list (optional).

    Returns:
        ``None``.

    """

    ax = pylab.gca()

    all_spines = ["top", "bottom", "right", "left", "polar"]

    for spine in all_spines:
        if spine in hidden:
            ax.spines[spine].set_visible(False)
        else:
            try:
                ax.spines[spine].set_visible(True)
                ax.spines[spine].set_linewidth(px2pt(0.75))
            except KeyError:
                pass

    return


def init_pylab(font_kwargs={}):
    """Initialise and clean up the look and feel of the plotting area.

    Returns:
        ``None``.

    """

    mpl.rc("lines", linewidth=px2pt(1))
    mpl.rc("xtick", **{"direction" : "out" })
    mpl.rc("ytick", **{"direction" : "out" })
    mpl.rc("legend", frameon=False, fontsize=font_kwargs["size"], numpoints=1)
    mpl.rc("font", **font_kwargs)

    pylab.tick_params(axis="x", which="both", top="off")
    pylab.tick_params(axis="y", which="both", right="off")

    init_spines()

    return


#
# Geometry
#

def calc_center_of_mass(atoms):
    """Compute the center of mass from a collection of atoms.

    Args:
        atoms: A list of Bio.PDB.Atom.Atom objects.

    Returns:
        The center of mass of ``atoms``.

    """

    coords = [a.get_coord() for a in atoms]
    weights = [a.mass for a in atoms]

    return numpy.average(coords, weights=weights, axis=0)


def calc_minvdw_distance(res_a, res_b):
    """Compute the minimum VDW distance between two residues, accounting for the
    VDW radii of each atom.

    Args:
        res_a: A ``Bio.PDB.Residue`` object.
        res_b: A ``Bio.PDB.Residue`` object.

    Returns:
        The minimum VDW distance between ``res_a`` and ``res_b``.

    """

    min_dist = None

    for a in res_a.get_iterator():
        for b in res_b.get_iterator():
            radii_a = VDW_RADII.get(a.get_id()[0], 0.0)
            radii_b = VDW_RADII.get(b.get_id()[0], 0.0)

            A = a.get_coord()
            B = b.get_coord()

            dist = numpy.linalg.norm(A - B) - radii_a - radii_b

            if (min_dist is None) or dist < min_dist:
                min_dist = dist

    return min_dist


def calc_cmass_distance(res_a, res_b, sidechain_only=False):
    """Compute the distance between the centres of mass of both residues.

    Args:
        res_a: A ``Bio.PDB.Residue`` object.
        res_b: A ``Bio.PDB.Residue`` object.
        sidechain_only: Set to True to consider only the sidechain atoms, False
            otherwise (optional).

    Returns:
        The distance between the centers of mass of ``res_a`` and ``res_b``.

    """

    if sidechain_only:
        atoms_a = get_sidechain_atoms(res_a)
        atoms_b = get_sidechain_atoms(res_b)
    else:
        atoms_a = res_a.get_list()
        atoms_b = res_b.get_list()

    A = calc_center_of_mass(atoms_a)
    B = calc_center_of_mass(atoms_b)

    return numpy.linalg.norm(A-B)


def calc_distance(res_a, res_b, measure="CA"):
    """Calculate the (L2) Euclidean distance between a pair of residues
    according to a given distance metric.

    Args:
        res_a: A ``Bio.PDB.Residue`` object.
        res_b: A ``Bio.PDB.Residue`` object.
        measure: The inter-residue distance measure (optional).

    Returns:
        The distance between ``res_a`` and ``res_b``.

    """

    if measure in ("CA", "CB"):
        A = get_atom_coord(res_a, measure)
        B = get_atom_coord(res_b, measure)
        dist = numpy.linalg.norm(A-B)
    elif measure == "cmass":
        dist = calc_cmass_distance(res_a, res_b)
    elif measure == "sccmass":
        dist = calc_cmass_distance(res_a, res_b, sidechain_only=True)
    elif measure == "minvdw":
        dist = calc_minvdw_distance(res_a, res_b)
    else:
        raise NotImplementedError

    return dist


def calc_dist_matrix(residues, measure="CA", dist_thresh=None,
        mask_thresh=None, asymmetric=False):
    """Calculate the distance matrix for a list of residues.

    Args:
        residues: A list of ``Bio.PDB.Residue`` objects.
        measure: The distance measure (optional).
        dist_thresh: (optional).
        mask_thresh: (optional).
        asymmetric: (optional).

    Returns:
        The distance matrix as a masked array.

    """

    mat = numpy.zeros((len(residues), len(residues)), dtype="float64")

    # after the distances are added to the upper-triangle, the nan values
    # indicate the lower matrix values, which are "empty", but can be used to
    # convey other information if needed.
    mat[:] = numpy.nan

    # Compute the upper-triangle of the underlying distance matrix.
    #
    # TODO:
    # - parallelise this over multiple processes + show benchmark results.
    # - use the lower-triangle to convey other information.
    pair_indices = combinations_with_replacement(range(len(residues)), 2)

    for i, j in pair_indices:
        res_a = residues[i]
        res_b = residues[j]
        dist = calc_distance(res_a, res_b, measure)
        mat[j,i] = dist

        if not asymmetric:
            mat[i,j] = dist

    # transpose i with j so the distances are contained only in the
    # upper-triangle.
    mat = mat.T
    mat = numpy.ma.masked_array(mat, numpy.isnan(mat))

    if dist_thresh is not None:
        mat = mat < dist_thresh

    if mask_thresh:
        mat = numpy.ma.masked_greater_equal(mat, mask_thresh)

    return mat



if __name__ == '__main__':
    opts = docopt(__doc__)

    # Distance threshold for contact maps.
    if opts["<dist>"]:
        opts["<dist>"] = float(opts["<dist>"])

    if opts["--mask-thresh"]:
        opts["--mask-thresh"] = float(opts["--mask-thresh"])

    if opts["--chains"]:
        chain_ids = opts["--chains"].upper().split(",")

        # Check that pdb chain ids are alphanumeric (see:
        # http://deposit.rcsb.org/adit/).
        if not numpy.all(list(map(str.isalnum, chain_ids))):
            sys.stderr.write()


    measure = opts["--measure"]

    residues = get_residues(opts["--pdb"], chain_ids=chain_ids)

    #
    # Generate the underlying 2D matrix for the selected plot.
    #
    mat = calc_dist_matrix(residues, measure=measure,
            dist_thresh=opts["<dist>"], mask_thresh=opts["--mask-thresh"],
            asymmetric=opts["--asymmetric"])

    if opts["--plaintext"]:
        if opts["cmap"]:
            fmt = "%d"
        else:
            fmt = "%.3f"

        numpy.savetxt(opts["--output"], mat.filled(0), fmt=fmt)
    else:
        font_kwargs = {
                "family" : opts["--font-family"],
                "size" : float(opts["--font-size"]) }

        init_pylab(font_kwargs)

        # hide all the spines i.e. no axes are drawn
        init_spines(hidden=["top", "bottom", "left", "right"])

        pylab.gcf().set_figwidth(float(opts["--width-inches"]))
        pylab.gcf().set_figheight(float(opts["--height-inches"]))

        pylab.xlim(0, len(residues))
        pylab.ylim(0, len(residues))

        pylab.xlabel(opts["--xlabel"])
        pylab.ylabel(opts["--ylabel"])

        ax, fig = pylab.gca(), pylab.gcf()

        if opts["--show-frame"]:
            init_spines(hidden=[])

        if opts["cmap"]:
            map_obj = pylab.pcolormesh(mat,
                    shading="flat", edgecolors="None", cmap=mpl.cm.Greys)
        elif opts["dmap"]:
            if opts["--greyscale"]:
                cmap = mpl.cm.Greys
            else:
                cmap = mpl.cm.jet_r

            map_obj = pylab.pcolormesh(mat, shading="flat",
                    edgecolors="None", cmap=cmap)

            if not opts["--no-colorbar"]:
                # draw the colour bar
                box = ax.get_position()
                pad, width = 0.02, 0.02
                cax = fig.add_axes([box.xmax + pad, box.ymin, width, box.height])
                cbar = pylab.colorbar(map_obj, drawedges=False, cax=cax)
                cbar.outline.set_visible(False)
                pylab.ylabel("Distance (angstroms)")
        else:
            raise NotImplementedError

        if opts["--title"] is not None:
            ax.set_title(opts["--title"], fontweight="bold")

        pylab.savefig(opts["--output"], bbox_inches="tight",
                dpi=int(opts["--dpi"]), transparent=opts["--transparent"])
