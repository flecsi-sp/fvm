#ifndef MUSCL_STATE_HH
#define MUSCL_STATE_HH

#include "types.hh"

namespace muscl {

inline mesh::slot m;

// Primitives.
inline const field<double>::definition<mesh, mesh::cells> u;
inline const field<double>::definition<mesh, mesh::cells> v;
inline const field<double>::definition<mesh, mesh::cells> w;
inline const field<double>::definition<mesh, mesh::cells> p;

// Intermediate conserved quantities.
inline const field<double>::definition<mesh, mesh::cells> q;
inline const field<double>::definition<mesh, mesh::cells> qu;
inline const field<double>::definition<mesh, mesh::cells> qv;
inline const field<double>::definition<mesh, mesh::cells> qw;
inline const field<double>::definition<mesh, mesh::cells> qE;

// Slopes.
inline const field<double>::definition<mesh, mesh::cells> dr_ds;
inline const field<double>::definition<mesh, mesh::cells> du_ds;
inline const field<double>::definition<mesh, mesh::cells> dv_ds;
inline const field<double>::definition<mesh, mesh::cells> dw_ds;
inline const field<double>::definition<mesh, mesh::cells> dp_ds;

// Faces.
inline const field<double>::definition<mesh, mesh::cells> rTail;
inline const field<double>::definition<mesh, mesh::cells> ruTail;
inline const field<double>::definition<mesh, mesh::cells> rvTail;
inline const field<double>::definition<mesh, mesh::cells> rETail;
inline const field<double>::definition<mesh, mesh::cells> uTail;
inline const field<double>::definition<mesh, mesh::cells> vTail;
inline const field<double>::definition<mesh, mesh::cells> pTail;

inline const field<double>::definition<mesh, mesh::cells> rHead;
inline const field<double>::definition<mesh, mesh::cells> ruHead;
inline const field<double>::definition<mesh, mesh::cells> rvHead;
inline const field<double>::definition<mesh, mesh::cells> rEHead;
inline const field<double>::definition<mesh, mesh::cells> uHead;
inline const field<double>::definition<mesh, mesh::cells> vHead;
inline const field<double>::definition<mesh, mesh::cells> pHead;

// Riemann fluxes.
inline const field<double>::definition<mesh, mesh::cells> rF;
inline const field<double>::definition<mesh, mesh::cells> ruF;
inline const field<double>::definition<mesh, mesh::cells> rvF;
inline const field<double>::definition<mesh, mesh::cells> rEF;

} // namespace muscl

#endif // MUSCL_STATE_HH
