// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef CRONOS__FFBLOCK_HPP
#define CRONOS__FFBLOCK_HPP

#include <algorithm>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "ffdep.hpp"
#include "ffunc.hpp"

// #define CRONOS__FFBLOCK_DEBUG
// #define CRONOS__FFBLOCK_DISABLE_MC21A

#ifdef CRONOS__WITH_HSL
// HSL MC21/MC13 (block triangular form) and MC33 (bordered block triangular
// form) -- licensed
extern "C" void mc13d_(const int*, const int*, const int*, const int*,
                       const int*, int*, int*, int*, int*);
extern "C" void mc21a_(const int*, const int*, const int*, const int*,
                       const int*, int*, int*, int*);
extern "C" void mc33ad_(const int*, const int*, const int*, int*, const int*,
                        double*, int*, int*, int*, int*, int*, int*, int*, int*,
                        int*);
#else
// SuiteSparse BTF (LGPL-2.1-or-later), distributed with KLU -- link with -lbtf
#include "btf.h"
#endif

namespace mc
{

//! @brief C++ class for block decompositions of equation systems in a DAG
////////////////////////////////////////////////////////////////////////
//! mc::FFBLOCK is a C++ class computing block decompositions of the
//! sparsity pattern of a system of equations defined in a DAG of MC++,
//! in CRONOS (so that MC++ has no HSL dependency):
//! - block_triangular: block lower-triangular form of a square system,
//!   with HSL MC21+MC13 if CRONOS__WITH_HSL is defined, else with SuiteSparse
//!   BTF (btf_maxtrans+btf_strongcomp, the same algorithms); the blocks
//!   are the strongly connected components in both cases
//! - bordered_triangular: bordered lower-triangular form of a possibly
//!   non-square system, with HSL MC33 (P5 algorithm, full diagonal blocks)
//!   if CRONOS__WITH_HSL is defined, else with a maximum transversal from BTF
//!   followed by a greedy tearing of the cycles (1x1 diagonal blocks); the
//!   two need not give the same border
//! All indices in the outputs are 1-based, as in HSL.
////////////////////////////////////////////////////////////////////////
class FFBLOCK
////////////////////////////////////////////////////////////////////////
{
 public:
  //! @brief Default Constructor
  FFBLOCK(FFGraph* dag = nullptr) : _dag(dag) {}

  //! @brief Destructor
  virtual ~FFBLOCK() {}

  //! @brief Retrieve pointer to DAG
  FFGraph*
  dag() const
  {
    return _dag;
  }

  //! @brief Set DAG environment
  void
  set(FFGraph* dag)
  {
    _dag = dag;
  }

  //! @brief Library used for the decompositions: "HSL" or "BTF"
  static char const*
  library()
  {
#ifdef CRONOS__WITH_HSL
    return "HSL";
#else
    return "BTF";
#endif
  }

  //! @brief Block lower-triangular rearrangement of the square system
  //! <a>pDep</a> in the variables <a>pIndep</a> -- equation IPERM[IOR[i]-1] and
  //! variable IOR[i] are in position i of the permuted system, block k
  //! (k=1..NB) starts at position IB[k-1], and block k depends on blocks 1..k
  //! only, as in HSL MC21 and MC13; false if the system is structurally
  //! singular
  bool block_triangular(unsigned const nDep, FFVar const* pDep,
                        FFVar const* pIndep, int* IPERM, int* IOR, int* IB,
                        int& NB, bool const disp = false,
                        std::ostream& os = std::cout);

  //! @brief Same as block_triangular for the sparsity pattern <a>rows</a> of a
  //! square system of size <a>n</a>: equation i depends on the variables
  //! rows[i] (0-based)
  static bool block_triangular(unsigned const n,
                               std::vector<std::vector<int>> const& rows,
                               int* IPERM, int* IOR, int* IB, int& NB);

  //! @brief Bordered lower-triangular rearrangement of the <a>nDep</a>
  //! equations <a>pDep</a> in the <a>nIndep</a> variables <a>pIndep</a> --
  //! equation IP[i] and variable IQ[j] are in row i and column j of the
  //! permuted system; rows NROW+1.. and columns NCOL+1.. are null; the leading
  //! block (rows 1..NROW-NBROW, columns 1..NCOL-NBCOL) is block
  //! lower-triangular, the border is rows NROW-NBROW+1..NROW and columns
  //! NCOL-NBCOL+1..NCOL; false on an input error
  bool bordered_triangular(unsigned const nDep, FFVar const* pDep,
                           unsigned const nIndep, FFVar const* pIndep, int* IP,
                           int* IQ, int& NROW, int& NCOL, int& NBROW,
                           int& NBCOL, bool const disp = false,
                           std::ostream& os = std::cout);

  //! @brief Same as bordered_triangular for the sparsity pattern <a>rows</a> of
  //! an <a>m</a>-by-<a>n</a> system: equation i depends on the variables
  //! rows[i] (0-based)
  static bool bordered_triangular(unsigned const m, unsigned const n,
                                  std::vector<std::vector<int>> const& rows,
                                  int* IP, int* IQ, int& NROW, int& NCOL,
                                  int& NBROW, int& NBCOL);

  //! @brief FFBLOCK Exceptions
  class Exceptions
  {
   public:
    //! @brief Enumeration type for exception handling
    enum TYPE
    {
      DAG   = 1,   //!< DAG not set
      UNDEF = -33  //!< Undocumented error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions(TYPE ierr = UNDEF) : _ierr(ierr) {}
    //! @brief Inline function returning the error flag
    int
    ierr()
    {
      return _ierr;
    }
    //! @brief Error description
    std::string
    what()
    {
      switch (_ierr)
      {
        case DAG:
          return "mc::FFBLOCK::Exceptions\t DAG not set";
        case UNDEF:
        default:
          return "mc::FFBLOCK::Exceptions\t Undocumented error";
      }
    }

   private:
    TYPE _ierr;
  };

 protected:
  //! @brief pointer to underlying dag
  FFGraph* _dag;

  //! @brief Sparsity pattern of the equations <a>pDep</a> in the variables
  //! <a>pIndep</a>, the other variables in the DAG treated as constants
  std::vector<std::vector<int>> _pattern(unsigned const nDep, FFVar const* pDep,
                                         unsigned const nIndep,
                                         FFVar const* pIndep) const;

 private:
  //! @brief Private methods to block default compiler methods
  FFBLOCK(FFBLOCK const&)            = delete;
  FFBLOCK& operator=(FFBLOCK const&) = delete;
};

////////////////////////////////// FFBLOCK
//////////////////////////////////////////

inline std::vector<std::vector<int>>
FFBLOCK::_pattern(unsigned const nDep, FFVar const* pDep, unsigned const nIndep,
                  FFVar const* pIndep) const
{
  if (!_dag) throw Exceptions(Exceptions::DAG);

  // Variables, with the other DAG variables as constants
  std::vector<FFVar> pVar(pIndep, pIndep + nIndep);
  std::vector<FFDep> vVar(nIndep);
  for (unsigned i = 0; i < nIndep; i++) vVar[i].indep(i);
  auto&& sgDep = _dag->subgraph(nDep, pDep);
  for (auto const& op : sgDep.l_op)
  {
    // Operation not a variable
    if (op->type != FFOp::VAR) continue;
    FFVar var    = *op->varout.front();
    bool isParam = true;
    // Operation is a current independent
    for (unsigned i = 0; isParam && i < nIndep; i++)
      if (var.id() == pIndep[i].id()) isParam = false;
    if (!isParam) continue;
    // Add dummy variable for parameter
    pVar.push_back(var);
    vVar.push_back(FFDep());
  }
  std::vector<FFDep> wkDep(sgDep.len_tap);
  std::vector<FFDep> vDep(nDep);
  _dag->eval(sgDep, wkDep, nDep, pDep, vDep.data(), pVar.size(), pVar.data(),
             vVar.data());

  std::vector<std::vector<int>> rows(nDep);
  for (unsigned i = 0; i < nDep; i++)
    for (auto const& [j, typ] : vDep[i].dep()) rows[i].push_back(j);
  return rows;
}

inline bool
FFBLOCK::block_triangular(unsigned const n,
                          std::vector<std::vector<int>> const& rows, int* IPERM,
                          int* IOR, int* IB, int& NB)
{
  NB = 0;
  if (!n) return true;

#ifdef CRONOS__WITH_HSL
  // Row-wise pattern (!!Fortran style indices!!)
  std::vector<int> IP(n), LENR(n), ICN;
  for (unsigned i = 0; i < n; i++)
  {
    IP[i]   = ICN.size() + 1;
    LENR[i] = rows[i].size();
    for (int j : rows[i]) ICN.push_back(j + 1);
  }

  // Make a row permutation to remove nonzeros on diagonal: MC21A
  int N = n, LICN = ICN.size();
  std::vector<int> IW(4 * n);
#ifdef CRONOS__FFBLOCK_DISABLE_MC21A
  for (unsigned i = 0; i < n; i++) IPERM[i] = i + 1;
#else
  int NUMNZ;
  mc21a_(&N, ICN.data(), &LICN, IP.data(), LENR.data(), IPERM, &NUMNZ,
         IW.data());
  if (NUMNZ < N) return false;  // structural singularity

  // Permute order of equation system using IPERM
  ICN.clear();
  for (unsigned i = 0; i < n; i++)
  {
    IP[i]   = ICN.size() + 1;
    LENR[i] = rows[IPERM[i] - 1].size();
    for (int j : rows[IPERM[i] - 1]) ICN.push_back(j + 1);
  }
#endif

  // Make a block lower-triangular decomposition: MC13D
  mc13d_(&N, ICN.data(), &LICN, IP.data(), LENR.data(), IOR, IB, &NB,
         IW.data());

#else
  // Compressed-column form: column j (variable) holds the equations depending
  // on it
  std::vector<int32_t> Ap(n + 1, 0), Ai;
  for (auto const& row : rows)
    for (int j : row) Ap[j + 1]++;
  for (unsigned j = 0; j < n; j++) Ap[j + 1] += Ap[j];
  Ai.resize(Ap[n] ? Ap[n] : 1);
  std::vector<int32_t> pos(Ap.begin(), Ap.end() - 1);
  for (unsigned i = 0; i < n; i++)
    for (int j : rows[i]) Ai[pos[j]++] = i;

  // Maximum transversal and strongly connected components: P*A*Q block
  // upper-triangular
  std::vector<int32_t> P(n), Q(n), R(n + 1), Work(5 * n);
  double work        = 0.;
  int32_t nmatch     = 0;
  int32_t const nblk = btf_order(n, Ap.data(), Ai.data(), 0., &work, P.data(),
                                 Q.data(), R.data(), &nmatch, Work.data());
  if (nmatch < (int32_t)n) return false;  // structural singularity

  // Reverse to the block lower-triangular form of MC13
  for (unsigned k = 0; k < n; k++)
  {
    IOR[n - 1 - k] = Q[k] + 1;  // variable in position n-1-k
    IPERM[Q[k]]    = P[k] + 1;  // equation matched with that variable
  }
  for (int32_t b = 0; b < nblk; b++) IB[nblk - 1 - b] = n - R[b + 1] + 1;
  NB = nblk;
#endif

#ifdef CRONOS__FFBLOCK_DEBUG
  std::cout << "Number of blocks in permuted matrix: " << NB << std::endl;
  std::cout << "Row/Column reordering: ";
  for (unsigned i = 0; i < n; i++) std::cout << " " << IOR[i];
  std::cout << std::endl;
#endif
  return true;
}

inline bool
FFBLOCK::block_triangular(unsigned const nDep, FFVar const* pDep,
                          FFVar const* pIndep, int* IPERM, int* IOR, int* IB,
                          int& NB, bool const disp, std::ostream& os)
{
  auto const rows = _pattern(nDep, pDep, nDep, pIndep);
  if (!block_triangular(nDep, rows, IPERM, IOR, IB, NB)) return false;

  // Display permuted system structure
  if (disp)
  {
    os << std::endl
       << "Number of Blocks: " << NB << " (" << library() << ")" << std::endl
       << "Lower-triangular block structure:" << std::endl
       << std::right << "     ";
    for (unsigned j = 0; j < nDep; j++)
      os << " " << std::setw(4) << pIndep[IOR[j] - 1];
    os << std::endl;
    for (unsigned i = 0; i < nDep; i++)
    {
      auto const& row = rows[IPERM[IOR[i] - 1] - 1];
      os << std::setw(4) << pDep[IPERM[IOR[i] - 1] - 1] << " ";
      for (unsigned j = 0; j < nDep; j++)
        os << std::setw(3) << " "
           << (std::find(row.begin(), row.end(), IOR[j] - 1) != row.end()
                   ? "X "
                   : "  ");
      os << std::endl;
    }
    os << std::endl;
  }
  return true;
}

inline bool
FFBLOCK::bordered_triangular(unsigned const m, unsigned const n,
                             std::vector<std::vector<int>> const& rows, int* IP,
                             int* IQ, int& NROW, int& NCOL, int& NBROW,
                             int& NBCOL)
{
  NROW = NCOL = NBROW = NBCOL = 0;
  if (!m || !n) return false;

#ifdef CRONOS__WITH_HSL
  // MC33 is called on the non-null rows and columns only: given null rows or
  // columns, it can insert the null row M (or column N) among the non-null
  // ones, although its specification places them last -- the null rows and
  // columns are appended afterwards
  std::vector<int> rowid, colid,
      colnew(n, -1);  // compressed -> original, original -> compressed
  for (unsigned i = 0; i < m; i++)
    if (!rows[i].empty()) rowid.push_back(i);
  for (unsigned i = 0; i < m; i++)
    for (int j : rows[i]) colnew[j] = 0;
  for (unsigned j = 0; j < n; j++)
    if (!colnew[j])
    {
      colnew[j] = colid.size();
      colid.push_back(j);
    }
  if (rowid.empty()) return false;  // no entries, as MC33

  // Coordinate form (!!Fortran style indices!!)
  std::vector<int> IRN, JCN;
  std::vector<double> A;
  for (unsigned k = 0; k < rowid.size(); k++)
    for (int j : rows[rowid[k]])
    {
      IRN.push_back(k + 1);
      JCN.push_back(colnew[j] + 1);
      A.push_back(IRN.size());
    }

  // Make a bordered-block lower-triangular decomposition: MC33A, P5 algorithm
  const int ITYPE = 5;
  const int M = rowid.size(), N = colid.size(), NZI = IRN.size();
  int NZO, IERR, IFLAG[3];
  std::vector<int> IPc(M), IQc(N), IPROF(N), IW(M + N), IW1(9 * N + 3 * M);
  mc33ad_(&M, &N, &NZI, &NZO, &ITYPE, A.data(), IRN.data(), JCN.data(),
          IPc.data(), IQc.data(), IPROF.data(), IFLAG, IW.data(), IW1.data(),
          &IERR);
  if (IERR || IFLAG[0] != M || IFLAG[1] != N) return false;
  NROW  = M;                 // last non-null row
  NCOL  = N;                 // last non-null column
  NBROW = NBCOL = IFLAG[2];  // dimension of the border

  // Permutations of the original system, null rows and columns last
  std::vector<char> rowused(m, 0), colused(n, 0);
  for (int k = 0; k < M; k++)
  {
    IP[k]                      = rowid[IPc[k] - 1] + 1;
    rowused[rowid[IPc[k] - 1]] = 1;
  }
  for (int k = 0; k < N; k++)
  {
    IQ[k]                      = colid[IQc[k] - 1] + 1;
    colused[colid[IQc[k] - 1]] = 1;
  }
  for (unsigned i = 0, k = M; i < m; i++)
    if (!rowused[i]) IP[k++] = i + 1;
  for (unsigned j = 0, k = N; j < n; j++)
    if (!colused[j]) IQ[k++] = j + 1;

#else
  // Maximum transversal: Match[i] = j if equation i is matched with variable j
  std::vector<int32_t> Ap(n + 1, 0), Ai;
  for (auto const& row : rows)
    for (int j : row) Ap[j + 1]++;
  for (unsigned j = 0; j < n; j++) Ap[j + 1] += Ap[j];
  if (!Ap[n]) return false;  // no entries, as MC33
  Ai.resize(Ap[n]);
  std::vector<int32_t> pos(Ap.begin(), Ap.end() - 1);
  for (unsigned i = 0; i < m; i++)
    for (int j : rows[i]) Ai[pos[j]++] = i;
  std::vector<int32_t> Match(m), Work(5 * n);
  double work = 0.;
  btf_maxtrans(m, n, Ap.data(), Ai.data(), 0., &work, Match.data(),
               Work.data());

  // Matched pairs k = (equation, variable) and the dependencies between them:
  // pair k depends on pair l if the equation of k involves the variable of l
  std::vector<int> pairrow, paircol, colpair(n, -1);
  for (unsigned i = 0; i < m; i++)
  {
    if (Match[i] < 0) continue;
    colpair[Match[i]] = pairrow.size();
    pairrow.push_back(i);
    paircol.push_back(Match[i]);
  }
  unsigned const np = pairrow.size();
  std::vector<std::vector<int>> dep(np), user(np);  // dependencies / dependents
  for (unsigned k = 0; k < np; k++)
    for (int j : rows[pairrow[k]])
    {
      int const l = colpair[j];
      if (l < 0 || l == (int)k) continue;
      dep[k].push_back(l);
      user[l].push_back(k);
    }
  for (unsigned k = 0; k < np; k++)
  {
    std::sort(dep[k].begin(), dep[k].end());
    dep[k].erase(std::unique(dep[k].begin(), dep[k].end()), dep[k].end());
    std::sort(user[k].begin(), user[k].end());
    user[k].erase(std::unique(user[k].begin(), user[k].end()), user[k].end());
  }

  // Greedy tearing: a pair with no remaining dependency goes next at the front,
  // a pair with no remaining dependent goes next at the back; otherwise the
  // pair with most cycles through it, estimated by
  // (#dependencies)x(#dependents), is torn -- its variable and equation go to
  // the border. The pairs kept are then in lower-triangular order.
  std::vector<int> ndep(np), nuser(np), front, back, torn;
  std::vector<char> alive(np, 1);
  for (unsigned k = 0; k < np; k++)
  {
    ndep[k]  = dep[k].size();
    nuser[k] = user[k].size();
  }
  auto remove = [&](int k)
  {
    alive[k] = 0;
    for (int l : dep[k])
      if (alive[l]) nuser[l]--;
    for (int l : user[k])
      if (alive[l]) ndep[l]--;
  };
  for (unsigned nalive = np; nalive;)
  {
    int kf = -1, kb = -1, kt = -1;
    long best = -1;
    for (unsigned k = 0; k < np && kf < 0; k++)
    {
      if (!alive[k]) continue;
      if (!ndep[k])
        kf = k;
      else if (!nuser[k] && kb < 0)
        kb = k;
      else if (kb < 0)
      {
        long const score = (long)ndep[k] * nuser[k];
        if (score > best)
        {
          best = score;
          kt   = k;
        }
      }
    }
    if (kf >= 0)
    {
      front.push_back(kf);
      remove(kf);
    }
    else if (kb >= 0)
    {
      back.push_back(kb);
      remove(kb);
    }
    else
    {
      torn.push_back(kt);
      remove(kt);
    }
    --nalive;
  }

  // Assemble permutations: kept pairs, unmatched non-null equations and
  // variables in the leading block or the border, null equations and variables
  // last
  std::vector<char> rowused(m, 0), colused(n, 0), colnull(n, 1);
  for (unsigned i = 0; i < m; i++)
    for (int j : rows[i]) colnull[j] = 0;
  unsigned ir = 0, ic = 0;
  auto putrow = [&](int i)
  {
    IP[ir++]   = i + 1;
    rowused[i] = 1;
  };
  auto putcol = [&](int j)
  {
    IQ[ic++]   = j + 1;
    colused[j] = 1;
  };
  for (int k : front)
  {
    putrow(pairrow[k]);
    putcol(paircol[k]);
  }
  for (auto it = back.rbegin(); it != back.rend(); ++it)
  {
    putrow(pairrow[*it]);
    putcol(paircol[*it]);
  }
  // unmatched equations with entries in leading variables only stay below the
  // leading block
  unsigned const nlead = ic;
  for (unsigned i = 0; i < m; i++)
  {
    if (rowused[i] || rows[i].empty() || Match[i] >= 0) continue;
    bool lead = true;
    for (int j : rows[i])
      if (!colused[j]) lead = false;
    if (lead) putrow(i);
  }
  unsigned const nleadrow = ir;
  for (int k : torn)
  {
    putrow(pairrow[k]);
    putcol(paircol[k]);
  }
  for (unsigned i = 0; i < m; i++)
    if (!rowused[i] && !rows[i].empty()) putrow(i);
  for (unsigned j = 0; j < n; j++)
    if (!colused[j] && !colnull[j]) putcol(j);
  NROW  = ir;
  NCOL  = ic;
  NBROW = NROW - nleadrow;
  NBCOL = NCOL - nlead;
  for (unsigned i = 0; i < m; i++)
    if (!rowused[i]) putrow(i);
  for (unsigned j = 0; j < n; j++)
    if (!colused[j]) putcol(j);
#endif

#ifdef CRONOS__FFBLOCK_DEBUG
  std::cout << "Border: " << NBROW << " rows, " << NBCOL << " columns"
            << std::endl;
#endif
  return true;
}

inline bool
FFBLOCK::bordered_triangular(unsigned const nDep, FFVar const* pDep,
                             unsigned const nIndep, FFVar const* pIndep,
                             int* IP, int* IQ, int& NROW, int& NCOL, int& NBROW,
                             int& NBCOL, bool const disp, std::ostream& os)
{
  auto const rows = _pattern(nDep, pDep, nIndep, pIndep);
  if (!bordered_triangular(nDep, nIndep, rows, IP, IQ, NROW, NCOL, NBROW,
                           NBCOL))
    return false;

  // Display permuted system structure
  if (disp)
  {
    os << std::endl
       << "Bordered lower-triangular structure (" << library()
       << "), border of " << NBROW << " rows and " << NBCOL
       << " columns:" << std::endl
       << std::right << "    ";
    for (unsigned j = 0; j < nIndep; j++)
      os << ((int)j == NCOL - NBCOL ? " | " : " ") << std::setw(4)
         << pIndep[IQ[j] - 1];
    os << std::endl;
    for (unsigned i = 0; i < nDep; i++)
    {
      if ((int)i == NROW - NBROW)
        os << "     " << std::string(5 * nIndep + 2, '-') << std::endl;
      auto const& row = rows[IP[i] - 1];
      os << std::setw(4) << pDep[IP[i] - 1];
      for (unsigned j = 0; j < nIndep; j++)
        os << ((int)j == NCOL - NBCOL ? "  | " : "   ")
           << (std::find(row.begin(), row.end(), IQ[j] - 1) != row.end()
                   ? "X "
                   : "  ");
      os << std::endl;
    }
    os << std::endl;
  }
  return true;
}

}  // namespace mc

#endif
