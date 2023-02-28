
============================

(Mostly) Implementation-agnostic SIMD trigonometry functions

:::   Exp    :::
:::   Log    :::
:::   Sin    :::
:::   Cos    :::
:::   Tan    :::
:::   Cot    :::
:::   Atan   :::
:::   Atan2  :::

Uses SSE2 semantics. No special optimizations made for
particular SIMD instruction sets. This is intended to be
"good enough" for common use, like say, in a video game ;)

Implement "simd_trig_skeleton.h" using some set of
SIMD instructions, like SSE2 or NEON, or even a fallback.
Import your implementation into this file.

With a little effort, could probably be made to work with
register widths other than 128bits (probably just needs some
changes to the various constant array sizes)

SSE1 and MMX special casing has been removed entirely.
It's 2023 here, and that's old news, sorry.

Based on "MathFun" and "MathFunExt", see the respective
commentary below. I feel bad even putting my name on this,
as the real work was done by the previous authors,
but am doing so just to avoid licensing problems.

============================


License: zlib
See simd_trig.h for more details