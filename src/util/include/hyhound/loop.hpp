#pragma once

#include <hyhound/config.hpp>
#include <hyhound/openmp.h>

namespace hyhound {

enum class LoopDir {
    Forward,
    Backward,
};

[[gnu::always_inline]] inline void
foreach_chunked(index_t i_begin, index_t i_end, auto chunk_size,
                auto func_chunk, auto func_rem,
                LoopDir dir = LoopDir::Forward) {
    if (dir == LoopDir::Forward) {
        index_t i;
        for (i = i_begin; i + chunk_size <= i_end; i += chunk_size)
            func_chunk(i);
        index_t rem_i = i_end - i;
        if (rem_i > 0)
            func_rem(i, rem_i);
    } else {
        index_t rem_i = (i_end - i_begin) % chunk_size;
        index_t i     = i_end - rem_i;
        if (rem_i > 0)
            func_rem(i, rem_i);
        for (i -= chunk_size; i >= i_begin; i -= chunk_size)
            func_chunk(i);
    }
}

[[gnu::always_inline]] inline void
foreach_chunked_merged(index_t i_begin, index_t i_end, auto chunk_size,
                       auto func_chunk, LoopDir dir = LoopDir::Forward) {
    if (dir == LoopDir::Forward) {
        index_t i;
        for (i = i_begin; i + chunk_size <= i_end; i += chunk_size)
            func_chunk(i, chunk_size);
        index_t rem_i = i_end - i;
        if (rem_i > 0)
            func_chunk(i, rem_i);
    } else {
        index_t rem_i = (i_end - i_begin) % chunk_size;
        index_t i     = i_end - rem_i;
        if (rem_i > 0)
            func_chunk(i, rem_i);
        for (i -= chunk_size; i >= i_begin; i -= chunk_size)
            func_chunk(i, chunk_size);
    }
}

[[gnu::always_inline]] inline void
foreach_chunked_merged_parallel(index_t i_begin, index_t i_end, auto chunk_size,
                                auto func_chunk,
                                LoopDir dir = LoopDir::Forward) {
    const index_t rem_i = (i_end - i_begin) % chunk_size;
    if (dir == LoopDir::Forward) {
        HYHOUND_OMP(parallel) {
            HYHOUND_OMP(for nowait)
            for (index_t i = i_begin; i <= i_end - chunk_size; i += chunk_size)
                func_chunk(i, chunk_size);
            HYHOUND_OMP(single) {
                if (rem_i > 0)
                    func_chunk(i_end - rem_i, rem_i);
            }
        }
    } else {
        const index_t i_last = i_end - rem_i;
        HYHOUND_OMP(parallel) {
            HYHOUND_OMP(single nowait) {
                if (rem_i > 0)
                    func_chunk(i_last, rem_i);
            }
            HYHOUND_OMP(for)
            for (index_t i = i_last - chunk_size; i >= i_begin; i -= chunk_size)
                func_chunk(i, chunk_size);
        }
    }
}

} // namespace hyhound
