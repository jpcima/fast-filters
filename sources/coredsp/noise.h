#pragma once
#include <cstdint>

namespace coredsp {

class WhiteNoise {
  uint32_t state_;

 public:
  explicit WhiteNoise(uint32_t seed = 777);
  float tick();
};

inline WhiteNoise::WhiteNoise(uint32_t seed)
    : state_(seed) {}

inline float WhiteNoise::tick() {
  state_ = state_ * 1664525 + 1013904223;
  return int32_t(state_) * float(1.0 / INT32_MAX);
}

}  // namespace coredsp
