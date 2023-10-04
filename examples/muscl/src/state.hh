#ifndef MUSCL_STATE_HH
#define MUSCL_STATE_HH

#include "types.hh"

namespace muscl {

inline mesh::slot m;
inline index::slot ct; /* color topology */

// Time parameters.
inline const single<double>::definition<flecsi::topo::index> t;
inline const single<double>::definition<flecsi::topo::index> dt;
inline const single<velocity>::definition<flecsi::topo::index> lmax;

// Conserved quantities.
inline const field<double>::definition<mesh, mesh::cells> r;
inline const field<velocity>::definition<mesh, mesh::cells> ru;
inline const field<double>::definition<mesh, mesh::cells> rE;

// Primitives.
inline const field<velocity>::definition<mesh, mesh::cells> u;
inline const field<double>::definition<mesh, mesh::cells> p;

// Intermediate conserved quantities.
inline const field<double>::definition<mesh, mesh::cells> q;
inline const field<velocity>::definition<mesh, mesh::cells> qu;
inline const field<double>::definition<mesh, mesh::cells> qE;

// Slopes.
inline const field<double>::definition<mesh, mesh::cells> dr_ds;
inline const field<velocity>::definition<mesh, mesh::cells> du_ds;
inline const field<double>::definition<mesh, mesh::cells> dp_ds;

// Faces.
inline const field<double>::definition<mesh, mesh::cells> rTail;
inline const field<velocity>::definition<mesh, mesh::cells> ruTail;
inline const field<double>::definition<mesh, mesh::cells> rETail;
inline const field<velocity>::definition<mesh, mesh::cells> uTail;
inline const field<double>::definition<mesh, mesh::cells> pTail;

inline const field<double>::definition<mesh, mesh::cells> rHead;
inline const field<velocity>::definition<mesh, mesh::cells> ruHead;
inline const field<double>::definition<mesh, mesh::cells> rEHead;
inline const field<velocity>::definition<mesh, mesh::cells> uHead;
inline const field<double>::definition<mesh, mesh::cells> pHead;

// Riemann fluxes.
inline const field<double>::definition<mesh, mesh::cells> rF;
inline const field<velocity>::definition<mesh, mesh::cells> ruF;
inline const field<double>::definition<mesh, mesh::cells> rEF;

} // namespace muscl

#endif // MUSCL_STATE_HH
