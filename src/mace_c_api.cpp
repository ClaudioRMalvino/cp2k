/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__MACE)

#include <cassert>
#include <cstddef>
#include <span>
#include <string>

#include "mace.hpp" // from the symmetrix library (wcwitt/symmetrix)

#ifdef __cplusplus
extern "C" {
#endif

/*******************************************************************************
 * \brief Loads a MACE model from a symmetrix JSON file (symmetrix_extract_mace).
 * \author Claudio Malvino
 ******************************************************************************/
void mace_c_model_load(MACE **model_out, const char *filename) {
  assert(*model_out == NULL);
  *model_out = new MACE(std::string(filename));
}

/*******************************************************************************
 * \brief Releases a MACE model and all its resources.
 * \author Claudio Malvino
 ******************************************************************************/
void mace_c_model_release(MACE *model) { delete (model); }

/*******************************************************************************
 * \brief Returns the radial cutoff [model units] of a MACE model.
 * \author Claudio Malvino
 ******************************************************************************/
double mace_c_model_r_cut(const MACE *model) { return model->r_cut; }

/*******************************************************************************
 * \brief Returns the number of chemical elements the model was trained on.
 * \author Claudio Malvino
 ******************************************************************************/
int mace_c_model_num_elements(const MACE *model) { return model->num_elements; }

/*******************************************************************************
 * \brief Copies the model's atomic numbers (one per element) into dest.
 *        dest must hold num_elements ints.
 * \author Claudio Malvino
 ******************************************************************************/
void mace_c_model_atomic_numbers(const MACE *model, int *dest) {
  for (std::size_t i = 0; i < model->atomic_numbers.size(); i++) {
    dest[i] = model->atomic_numbers[i];
  }
}

/*******************************************************************************
 * \brief Copies the model's per-element reference (isolated-atom) energies
 *        into dest. dest must hold num_elements doubles.
 * \author Claudio Malvino
 ******************************************************************************/
void mace_c_model_atomic_energies(const MACE *model, double *dest) {
  for (std::size_t i = 0; i < model->atomic_energies.size(); i++) {
    dest[i] = model->atomic_energies[i];
  }
}

/*******************************************************************************
 * \brief Evaluates the MACE model over a local CSR neighbour graph.
 *
 *        Inputs (all caller-owned, raw pointers wrapped as std::span):
 *          num_nodes      number of local atoms (graph nodes)
 *          num_edges      total number of directed edges = sum(num_neigh)
 *          node_types     [num_nodes] MACE element index of each node
 *          num_neigh      [num_nodes] neighbour count of each node (CSR rows)
 *          neigh_indices  [num_edges] node index of each edge's neighbour
 *          neigh_types    [num_edges] MACE element index of each neighbour
 *          xyz            [3*num_edges] per-edge displacement vector
 *          r              [num_edges] per-edge distance
 *
 *        Outputs (caller-owned):
 *          node_energies_out  [num_nodes]   per-node energy
 *          node_forces_out    [3*num_edges] per-edge force (scatter to atoms
 *                                           is done on the Fortran side)
 * \author Claudio Malvino
 ******************************************************************************/
void mace_c_model_compute(MACE *model, const int num_nodes, const int num_edges,
                          const int *node_types, const int *num_neigh,
                          const int *neigh_indices, const int *neigh_types,
                          const double *xyz, const double *r,
                          double *node_energies_out, double *node_forces_out) {
  model->compute_node_energies_forces(
      num_nodes, std::span<const int>(node_types, num_nodes),
      std::span<const int>(num_neigh, num_nodes),
      std::span<const int>(neigh_indices, num_edges),
      std::span<const int>(neigh_types, num_edges),
      std::span<const double>(xyz, static_cast<std::size_t>(3) * num_edges),
      std::span<const double>(r, num_edges));

  for (int i = 0; i < num_nodes; i++) {
    node_energies_out[i] = model->node_energies[i];
  }
  for (int i = 0; i < 3 * num_edges; i++) {
    node_forces_out[i] = model->node_forces[i];
  }
}

#ifdef __cplusplus
}
#endif

#endif // defined(__MACE)

// EOF
