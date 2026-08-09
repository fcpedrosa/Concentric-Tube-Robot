#pragma once

// The Observer class that previously lived here has been removed.
//
// It stored ODE integration state via a functor (operator()) but was replaced
// by a capturing lambda directly inside CTR::ODESolver, which is simpler,
// avoids a dangling-reference risk, and removes the class from the public API.
//
// The state_type and bvp_type aliases are now in CTRTypes.hpp (ctr namespace).
