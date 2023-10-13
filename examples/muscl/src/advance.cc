#include "advance.hh"
#include "limiter.hh"
#include "state.hh"
#include "tasks/boundary.hh"
#include "tasks/hydro.hh"

#include <flecsi/execution.hh>

using namespace flecsi;
using namespace muscl;

void
action::advance(control_policy & cp) {
  // clang-format off
  execute<tasks::hydro::advance<mesh::x_axis, genminmod>>(m,
    r(m), ru(m), rE(m),
    u(m), p(m),
    q(m), qu(m), qE(m),
    dr_ds(m), du_ds(m), dp_ds(m),
    rTail(m), ruTail(m), rETail(m), uTail(m), pTail(m),
    rHead(m), ruHead(m), rEHead(m), uHead(m), pHead(m),
    rF(m), ruF(m), rEF(m),
    gamma(gt), cp.dt());
  // clang-format on

  execute<tasks::apply_boundaries>(m, bmap(gt), r(m), ru(m), rE(m));
  execute<tasks::hydro::update_primitives>(
    m, r(m), ru(m), rE(m), u(m), p(m), gamma(gt));

  // clang-format off
  execute<tasks::hydro::advance<mesh::y_axis, genminmod>>(m,
    r(m), ru(m), rE(m),
    u(m), p(m),
    q(m), qu(m), qE(m),
    dr_ds(m), du_ds(m), dp_ds(m),
    rTail(m), ruTail(m), rETail(m), uTail(m), pTail(m),
    rHead(m), ruHead(m), rEHead(m), uHead(m), pHead(m),
    rF(m), ruF(m), rEF(m),
    gamma(gt), cp.dt());
  // clang-format on

  execute<tasks::apply_boundaries>(m, bmap(gt), r(m), ru(m), rE(m));
  execute<tasks::hydro::update_primitives>(
    m, r(m), ru(m), rE(m), u(m), p(m), gamma(gt));

  // clang-format off
  execute<tasks::hydro::advance<mesh::z_axis, genminmod>>(m,
    r(m), ru(m), rE(m),
    u(m), p(m),
    q(m), qu(m), qE(m),
    dr_ds(m), du_ds(m), dp_ds(m),
    rTail(m), ruTail(m), rETail(m), uTail(m), pTail(m),
    rHead(m), ruHead(m), rEHead(m), uHead(m), pHead(m),
    rF(m), ruF(m), rEF(m),
    gamma(gt), cp.dt());
  // clang-format on

  execute<tasks::apply_boundaries>(m, bmap(gt), r(m), ru(m), rE(m));
  execute<tasks::hydro::update_primitives>(
    m, r(m), ru(m), rE(m), u(m), p(m), gamma(gt));

  execute<tasks::hydro::update_eigenvalues>(
    m, r(m), u(m), p(m), lmax(ct), gamma(gt));

  cp.dtmin() = reduce<tasks::hydro::update_dtmin, exec::fold::min>(m, lmax(ct));
} // action::advance
