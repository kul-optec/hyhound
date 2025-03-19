#pragma once

#define HYHOUND_CAT(a, b) HYHOUND_CAT_IMPL(a, b)
#define HYHOUND_CAT_IMPL(a, b) a##b
#define HYHOUND_STRINGIFY(s) HYHOUND_STRINGIFY_IMPL(s)
#define HYHOUND_STRINGIFY_IMPL(s) #s
