#include "advance.hh"
#include "limiter.hh"
#include "state.hh"
#include "tasks/advance.hh"

#include <flecsi/execution.hh>

using namespace flecsi;
using namespace muscl;

void
action::advance(control_policy & cp) {
  flog(info) << "advance action" << std::endl;

  const double dt{0};
  const double gamma{0};

  execute<tasks::advance<gminmod>>(m,
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
    dt,
    gamma);
} // action::advance
