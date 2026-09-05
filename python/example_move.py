"""
Example Python Move Script for ClassyMC

This script demonstrates how to create a custom Monte Carlo move
using Python. The move is called from Fortran via the forpy interface.

To use this move in your simulation, add the following to your input file:
    create moves 1
        python 1.0 example_move

Required Functions:
    - get_move_type(): Return the displacement type this move will generate
    - prologue(boxlist): Called once at the start of the simulation
    - generate_move(boxlist): Generate trial move, return displacement dict
    - compute_reverse_prob(boxlist, move_dict): Compute reverse probability
    - accept_move(boxlist, accepted, forward_prob, reverse_prob): Called after move decision

Move Types:
    - 'displacement': Single atom/molecule translation
    - 'deletion': Remove a molecule
    - 'addition': Add a molecule  
    - 'atomexchange': Exchange atom types
    - 'volchange': Isotropic volume change
    - 'orthovolchange': Anisotropic volume change
    - 'newstate_isomol': Same molecule count, many atom positions (full recompute)

The boxlist contains dictionaries with box information:
    - 'boxtype': string - type of box ('cube', 'ortho', etc.)
    - 'temperature': float - temperature of the box
    - 'pressure': float - pressure of the box  
    - 'volume': float - volume of the box
    - 'boxdimensions': numpy array - box dimensions
    - 'chemicalpotential': numpy array - chemical potentials
    - 'atomtype': numpy array - atom types
    - 'raw_atoms': numpy array - atomic coordinates (3 x nAtoms)
    - 'moleculecount': numpy array - number of molecules per type
    - 'energytable': numpy array - energy table

Probability Handling for Detailed Balance:
    - generate_move() can return 'forward_prob' in the dict (default 1.0)
    - compute_reverse_prob() returns the reverse generation probability
    - The acceptance criterion uses: Prob_ratio = reverse_prob / forward_prob
    - For symmetric moves (equal forward/reverse), both can be 1.0
    - For biased selection, set appropriate probabilities
"""

import numpy as np
import random

# Global state for the move (optional)
max_displacement = 0.5
accepted_count = 0
attempted_count = 0


def get_move_type():
    """
    Return the displacement type this move will generate.
    
    This function is called once during initialization to determine
    what type of perturbation array to allocate.
    
    Returns:
        str: One of:
            - 'displacement' - single atom/molecule translation
            - 'deletion' - remove a molecule
            - 'addition' - add a molecule
            - 'atomexchange' - exchange atom types
            - 'volchange' - isotropic volume change
            - 'orthovolchange' - anisotropic volume change
            - 'newstate_isomol' - same N, many atom positions
    """
    return 'displacement'


def prologue(boxlist):
    """
    Called once at the start of the simulation.
    Use this to initialize any state needed by the move.
    
    Args:
        boxlist: List of box dictionaries containing simulation state
    """
    global max_displacement
    print("Python Move: Initializing example move")
    print(f"  Number of boxes: {len(boxlist)}")
    for i, box in enumerate(boxlist):
        print(f"  Box {i}: {box['boxtype']}, T={box['temperature']}")
    

def generate_move(boxlist):
    """
    Generate a trial Monte Carlo move.
    
    This function is called for each move attempt. It should:
    1. Select which atom(s) to move
    2. Generate new position(s)
    3. Return a dictionary with the displacement information
    
    The dictionary format depends on the move type returned by get_move_type():
    
    For 'displacement':
        {'moltype': int, 'molindex': int, 'atomindex': int,
         'x_new': float, 'y_new': float, 'z_new': float}
    
    For 'deletion':
        {'moltype': int, 'molindex': int, 'atomindex': int}
    
    For 'addition':
        {'moltype': int, 'molindex': int, 'atomindex': int,
         'x_new': float, 'y_new': float, 'z_new': float}
    
    For 'atomexchange':
        {'oldatomindex': int, 'oldtype': int, 
         'newatomindex': int, 'newtype': int}
    
    For 'volchange':
        {'volnew': float, 'volold': float}
    
    For 'orthovolchange':
        {'volnew': float, 'volold': float,
         'xscale': float, 'yscale': float, 'zscale': float}

    For 'newstate_isomol':
        {'new_atoms': numpy array (3 x nMaxAtoms),
         'n_moved': int (optional)}

    Probability (either linear or log, not both):
        {'forward_prob': float}  - linear generation probability (default 1.0)
        {'log_forward_prob': float} - if present, compute_reverse_prob must
          return a log reverse probability. Used by force-bias MC to avoid
          underflow. See example_fbmc_move.py.

    Any type can also return:
        {'reject': True} - to immediately reject the move
    
    Args:
        boxlist: List of box dictionaries containing current simulation state
    
    Returns:
        dict: Dictionary with move information
    """
    global max_displacement, attempted_count
    attempted_count += 1
    
    # Get the first box (for single-box simulations)
    box = boxlist[0]
    atoms = box['raw_atoms']
    natoms = atoms.shape[1]
    
    if natoms == 0:
        # No atoms to move
        return {'reject': True}
    
    # Select a random atom
    atom_idx = random.randint(0, natoms - 1)
    
    # Get current position
    x_old = atoms[0, atom_idx]
    y_old = atoms[1, atom_idx]
    z_old = atoms[2, atom_idx]
    
    # Generate random displacement
    dx = max_displacement * (2.0 * random.random() - 1.0)
    dy = max_displacement * (2.0 * random.random() - 1.0)
    dz = max_displacement * (2.0 * random.random() - 1.0)
    
    # Calculate new position
    x_new = x_old + dx
    y_new = y_old + dy
    z_new = z_old + dz
    
    # Get molecule information for this atom
    # Note: atomtype array gives the type, but we need molecule info
    # For simplicity, we'll use placeholder values here
    # In a real implementation, you'd need to track this properly
    moltype = 1  # Molecule type (1-indexed)
    molindex = 1  # Global molecule index (1-indexed)
    atomindex = atom_idx + 1  # Atom index (1-indexed for Fortran)
    
    # For this simple uniform random selection, forward_prob = 1/natoms
    # Since the reverse selection would also be 1/natoms (symmetric),
    # we can either set both to 1.0 or set both to 1/natoms - the ratio is 1.0
    forward_prob = 1.0 / natoms
    
    return {
        'moltype': moltype,
        'molindex': molindex,
        'atomindex': atomindex,
        'x_new': x_new,
        'y_new': y_new,
        'z_new': z_new,
        'forward_prob': forward_prob
    }


def compute_reverse_prob(boxlist, move_dict):
    """
    Compute the reverse generation probability for detailed balance.
    
    This function is called after generate_move() returns, before the
    acceptance decision. It receives the same boxlist (still in the OLD
    configuration) and the move_dict returned by generate_move().
    
    The reverse probability is the probability of generating the reverse
    move (going from the trial state back to the current state).
    
    For symmetric moves like simple random displacement:
        - Forward: select atom i with prob 1/N, displace uniformly
        - Reverse: select atom i with prob 1/N, displace back
        - Since selection is uniform and displacement is symmetric,
          forward_prob = reverse_prob, so ratio = 1.0
    
    For biased moves (e.g., preferential selection based on energy):
        - Forward: select atom i with prob P_fwd(i)
        - Reverse: select atom i with prob P_rev(i) 
        - These may differ if selection depends on configuration

    For force-bias / smart MC (see example_fbmc_move.py):
        Return log T from generate_move as 'log_forward_prob' and return
        log T_reverse from this function. Classy then uses
            logProb = log_reverse - log_forward
        in the Metropolis criterion.
    
    Args:
        boxlist: List of box dictionaries (current/old configuration)
        move_dict: Dictionary returned by generate_move()
    
    Returns:
        float: Reverse generation probability
    """
    # For this symmetric uniform selection, reverse_prob = forward_prob
    box = boxlist[0]
    atoms = box['raw_atoms']
    natoms = atoms.shape[1]
    
    if natoms == 0:
        return 1.0
    
    # Uniform random selection: reverse probability is also 1/N
    reverse_prob = 1.0 / natoms
    return reverse_prob


def accept_move(boxlist, accepted, forward_prob, reverse_prob):
    """
    Called after the move acceptance decision.
    
    Use this to update internal state based on whether the move was accepted.
    This is useful for adaptive moves that adjust parameters based on
    acceptance rates, or for tracking probability statistics.
    
    Args:
        boxlist: List of box dictionaries (already updated if move was accepted)
        accepted: bool - True if the move was accepted, False otherwise
        forward_prob: float - The forward generation probability that was used
        reverse_prob: float - The reverse generation probability that was used
    
    Note:
        The acceptance criterion used: 
            acc = min(1, (reverse_prob/forward_prob) * exp(-beta * dE))
        where dE is the energy difference.
    """
    global accepted_count, max_displacement, attempted_count
    
    if accepted:
        accepted_count += 1
    
    # Optionally adjust max_displacement based on acceptance rate
    # (Simple adaptive scheme - adjust every 100 attempts)
    if attempted_count > 0 and attempted_count % 100 == 0:
        acceptance_rate = accepted_count / attempted_count
        if acceptance_rate > 0.5:
            max_displacement *= 1.05  # Increase displacement
        elif acceptance_rate < 0.3:
            max_displacement *= 0.95  # Decrease displacement
        prob_ratio = reverse_prob / forward_prob if forward_prob > 0 else 1.0
        print(f"Python Move: Acceptance rate = {acceptance_rate:.2%}, "
              f"max_displacement = {max_displacement:.4f}, "
              f"last prob_ratio = {prob_ratio:.4f}")
