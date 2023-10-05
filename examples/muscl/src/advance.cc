#include "advance.hh"
#include "limiter.hh"
#include "state.hh"
#include "tasks/advance.hh"

#include <flecsi/execution.hh>

using namespace flecsi;
using namespace muscl;

void
action::advance(control_policy & cp) {
  execute<tasks::advance<genminmod>>(m,
    r(m),
    ru(m),
    rE(m),
    u(m),
    p(m),
    q(m),
    qu(m),
    qE(m),
    dr_ds(m),
    du_ds(m),
    dp_ds(m),
    rTail(m),
    ruTail(m),
    rETail(m),
    uTail(m),
    pTail(m),
    rHead(m),
    ruHead(m),
    rEHead(m),
    uHead(m),
    pHead(m),
    rF(m),
    ruF(m),
    rEF(m),
    gamma(gt),
    cp.dt());
  cp.dtmin() = reduce<tasks::dtmin, exec::fold::min>(m, lmax(ct));
} // action::advance
