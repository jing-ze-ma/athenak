//========================================================================================
// AthenaXXX astrophysical plasma code
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file two_stream_rt_chain_nn.cpp
//! \brief tusplit-1102: explicit instantiations of TsrtCkChain (and TsrtGreyRt) for part
//! of the
//! tag sets picket_fence_two_stream_RT_pass dispatches (generated; the full list is
//! the union of the two_stream_rt_chain_*.cpp files).  Tags: (SPH, BSP, CCH, FRM,
//! FOP, JAC); the radial tier NN is 0 for every set.

#include <tuple>
#include <type_traits>

#include "utils/two_stream_rt_chain.hpp"

namespace two_stream_rt {

template void TsrtCkChain(TsrtCkChainCtx &, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>,
    std::integral_constant<int, 2>, std::integral_constant<int, 0>,
    std::bool_constant<true>, std::bool_constant<false>);
template void TsrtCkChain(TsrtCkChainCtx &, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>,
    std::integral_constant<int, 1>, std::integral_constant<int, 0>,
    std::bool_constant<true>, std::bool_constant<false>);
template void TsrtCkChain(TsrtCkChainCtx &, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>,
    std::integral_constant<int, 0>, std::integral_constant<int, 0>,
    std::bool_constant<true>, std::bool_constant<false>);
template void TsrtCkChain(TsrtCkChainCtx &, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>,
    std::integral_constant<int, 2>, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>);
template void TsrtCkChain(TsrtCkChainCtx &, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>,
    std::integral_constant<int, 1>, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>);
template void TsrtCkChain(TsrtCkChainCtx &, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>,
    std::integral_constant<int, 0>, std::integral_constant<int, 0>,
    std::bool_constant<false>, std::bool_constant<false>);
template void TsrtGreyRt(TsrtGreyRtCtx &, std::integral_constant<int, 72>);
template void TsrtGreyRt(TsrtGreyRtCtx &, std::integral_constant<int, 136>);
template void TsrtGreyRt(TsrtGreyRtCtx &, std::integral_constant<int, 264>);
template void TsrtGreyRt(TsrtGreyRtCtx &, std::integral_constant<int, 520>);
template void TsrtGreyRt(TsrtGreyRtCtx &, std::integral_constant<int, 1032>);

}  // namespace two_stream_rt
