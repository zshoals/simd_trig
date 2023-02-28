

/* 	simd_trig.h

	(Mostly) Implementation-agnostic SIMD trigonometry functions

	Types:
	* float32x4 equivalent ---- simd_fq ("float quad")
	* int32x4 equivalent ------ simd_iq ("int quad")
	* uint32x4 equivalent ----- simd_uq ("unsigned int quad")

	Functions:
	* sin ---- simd_fq_sin
	* cos ---- simd_fq_cos
	* tan ---- simd_fq_tan
	* cot ---- simd_fq_cot
	* atan --- simd_fq_atan
	* atan2 -- simd_fq_atan2
	* exp ---- simd_fq_exp
	* log ---- simd_fq_log

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
*/

/* Copyright (C) 2023  zshoals
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)
*/







/* SIMD (SSE1+MMX or SSE2) implementation of sin, cos, exp and log
   Inspired by Intel Approximate Math library, and based on the
   corresponding algorithms of the cephes math library
   The default is to use the SSE1 version. If you define USE_SSE2 the
   the SSE2 intrinsics will be used in place of the MMX intrinsics. Do
   not expect any significant performance improvement with SSE2.
*/

/* Copyright (C) 2007  Julien Pommier
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)
*/



#pragma once

#include <stdint.h>

//==================================================
//IMPORTANT: Import your skeleton implementation here

#include "YOU_MUST_REPLACE_ME_WITH_YOUR_IMPLEMENTATION_OR_WHATEVER.h"

//==================================================

//Removed for portability, now we're unaligned everywhere
#define ALIGN16_BEG
#define ALIGN16_END

#define _PS_CONST(Name, Val)                                            \
  static const ALIGN16_BEG float _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PI32_CONST(Name, Val)                                            \
  static const ALIGN16_BEG int32_t _pi32_##Name[4] ALIGN16_END = { Val, Val, Val, Val }
#define _PS_CONST_TYPE(Name, Type, Val)                                 \
  static const ALIGN16_BEG Type _ps_##Name[4] ALIGN16_END = { Val, Val, Val, Val }

_PS_CONST(1  , 1.0f);
_PS_CONST(0p5, 0.5f);
/* the smallest non denormalized float number */
_PS_CONST_TYPE(min_norm_pos, int32_t, 0x00800000);
_PS_CONST_TYPE(mant_mask, int32_t, 0x7f800000);
_PS_CONST_TYPE(inv_mant_mask, int32_t, ~0x7f800000);

_PS_CONST_TYPE(sign_mask, int32_t, (int32_t)0x80000000);
_PS_CONST_TYPE(inv_sign_mask, int32_t, ~0x80000000);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);

_PS_CONST(cephes_SQRTHF, 0.707106781186547524f);
_PS_CONST(cephes_log_p0, 7.0376836292E-2f);
_PS_CONST(cephes_log_p1, - 1.1514610310E-1f);
_PS_CONST(cephes_log_p2, 1.1676998740E-1f);
_PS_CONST(cephes_log_p3, - 1.2420140846E-1f);
_PS_CONST(cephes_log_p4, + 1.4249322787E-1f);
_PS_CONST(cephes_log_p5, - 1.6668057665E-1f);
_PS_CONST(cephes_log_p6, + 2.0000714765E-1f);
_PS_CONST(cephes_log_p7, - 2.4999993993E-1f);
_PS_CONST(cephes_log_p8, + 3.3333331174E-1f);
_PS_CONST(cephes_log_q1, -2.12194440e-4f);
_PS_CONST(cephes_log_q2, 0.693359375f);

static inline simd_fq simd_fq_log(simd_fq x)
{
	simd_iq emm0;
	simd_fq one = fq_load_u(_ps_1);
	simd_fq invalid_mask = fq_cmple(x, fq_zeroes());

	x = fq_max(x, iq_as_fq(iq_load_u(_ps_min_norm_pos)));

	emm0 = uq_as_iq(uq_rshift_23(fq_as_uq(x)));

	x = fq_and(x, iq_as_fq(iq_load_u(_ps_inv_mant_mask)));
	x = fq_or(x, fq_load_u(_ps_0p5));

	emm0 = iq_sub(emm0, iq_load_u(_pi32_0x7f));
	simd_fq e = iq_convert_fq(emm0);

	e = fq_add(e, one);

	simd_fq mask = fq_cmplt(x, fq_load_u(_ps_cephes_SQRTHF));
	simd_fq tmp = fq_and(x, mask);
	x = fq_sub(x, one);
	e = fq_sub(e, fq_and(one, mask));
	x = fq_add(x, tmp);

	simd_fq z = fq_mul(x, x);

	simd_fq y = fq_load_u(_ps_cephes_log_p0);
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p1));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p2));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p3));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p4));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p5));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p6));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p7));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_log_p8));
	y = fq_mul(y, x);

	y = fq_mul(y, z);

	tmp = fq_mul(e, fq_load_u(_ps_cephes_log_q1));
	y = fq_add(y, tmp);

	tmp = fq_mul(z, fq_load_u(_ps_0p5));
	y = fq_sub(y, tmp);

	tmp = fq_mul(e, fq_load_u(_ps_cephes_log_q2));
	x = fq_add(x, y);
	x = fq_add(x, tmp);
	x = fq_or(x, invalid_mask);

	return x;
}







_PS_CONST(exp_hi,	88.3762626647949f);
_PS_CONST(exp_lo,	-88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341f);
_PS_CONST(cephes_exp_C1, 0.693359375f);
_PS_CONST(cephes_exp_C2, -2.12194440e-4f);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4f);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3f);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3f);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2f);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1f);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1f);


static inline simd_fq simd_fq_exp(simd_fq x)
{
	simd_fq tmp = fq_zeroes();
	simd_fq fx;
	simd_iq emm0;

	simd_fq one = fq_load_u(_ps_1);

	x = fq_min(x, fq_load_u(_ps_exp_hi));
	x = fq_max(x, fq_load_u(_ps_exp_lo));

	fx = fq_mul(x, fq_load_u(_ps_cephes_LOG2EF));
	fx = fq_add(fx, fq_load_u(_ps_0p5));

	emm0 = fq_truncate_iq(fx);
	tmp = iq_convert_fq(emm0);

	simd_fq mask = fq_cmpgt(tmp, fx);
	mask = fq_and(mask, one);
	fx = fq_sub(tmp, mask);

	tmp = fq_mul(fx, fq_load_u(_ps_cephes_exp_C1));
	simd_fq z = fq_mul(fx, fq_load_u(_ps_cephes_exp_C2));
	x = fq_sub(x, tmp);
	x = fq_sub(x, z);

	z = fq_mul(x, x);

	simd_fq y = fq_load_u(_ps_cephes_exp_p0);
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_exp_p1));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_exp_p2));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_exp_p3));
	y = fq_mul(y, x);
	y = fq_add(y, fq_load_u(_ps_cephes_exp_p4));
	y = fq_mul(y, z);
	y = fq_add(y, x);
	y = fq_add(y, one);

	emm0 = fq_truncate_iq(fx);
	emm0 = iq_add(emm0, iq_load_u(_pi32_0x7f));
	emm0 = uq_as_iq(uq_lshift_23(iq_as_uq(emm0)));

	simd_fq pow2n = iq_as_fq(emm0);
	y = fq_mul(y, pow2n);

	return y;
}



_PS_CONST(minus_cephes_DP1, -0.78515625f);
_PS_CONST(minus_cephes_DP2, -2.4187564849853515625e-4f);
_PS_CONST(minus_cephes_DP3, -3.77489497744594108e-8f);
_PS_CONST(sincof_p0, -1.9515295891E-4f);
_PS_CONST(sincof_p1,  8.3321608736E-3f);
_PS_CONST(sincof_p2, -1.6666654611E-1f);
_PS_CONST(coscof_p0,  2.443315711809948E-005f);
_PS_CONST(coscof_p1, -1.388731625493765E-003f);
_PS_CONST(coscof_p2,  4.166664568298827E-002f);
_PS_CONST(cephes_FOPI, 1.27323954473516f); // 4 / M_PI

static inline simd_fq simd_fq_sin(simd_fq x)
{
	simd_fq xmm1;
	simd_fq xmm2 = fq_zeroes();
	simd_fq xmm3;
	simd_fq sign_bit;
	simd_fq y;

	simd_iq emm0;
	simd_iq emm2;

	sign_bit = x;

	x = fq_and(x, iq_as_fq(iq_load_u(_ps_inv_sign_mask)));

	sign_bit = fq_and(sign_bit, iq_as_fq(iq_load_u(_ps_sign_mask)));

	y = fq_mul(x, fq_load_u(_ps_cephes_FOPI));

	emm2 = fq_truncate_iq(y);

	emm2 = iq_add(emm2, iq_load_u(_pi32_1));
	emm2 = iq_and(emm2, iq_load_u(_pi32_inv1));
	y = iq_convert_fq(emm2);

	emm0 = iq_and(emm2, iq_load_u(_pi32_4));
	emm0 = uq_as_iq(uq_lshift_29(iq_as_uq(emm0)));

	emm2 = iq_and(emm2, iq_load_u(_pi32_2));
	emm2 = uq_as_iq(uq_cmpeq(iq_as_uq(emm2), uq_zeroes()));

	simd_fq swap_sign_bit = iq_as_fq(emm0);
	simd_fq poly_mask = iq_as_fq(emm2);
	sign_bit = fq_xor(sign_bit, swap_sign_bit);

	xmm1 = fq_load_u(_ps_minus_cephes_DP1);
	xmm2 = fq_load_u(_ps_minus_cephes_DP2);
	xmm3 = fq_load_u(_ps_minus_cephes_DP3);
	xmm1 = fq_mul(y, xmm1);
	xmm2 = fq_mul(y, xmm2);
	xmm3 = fq_mul(y, xmm3);
	x = fq_add(x, xmm1);
	x = fq_add(x, xmm2);
	x = fq_add(x, xmm3);

	y = fq_load_u(_ps_coscof_p0);
	simd_fq z = fq_mul(x, x);

	y = fq_mul(y, z);
	y = fq_add(y, fq_load_u(_ps_coscof_p1));
	y = fq_mul(y, z);
	y = fq_add(y, fq_load_u(_ps_coscof_p2));
	//Yes, twice
	y = fq_mul(y, z);
	y = fq_mul(y, z);
	simd_fq tmp = fq_mul(z, fq_load_u(_ps_0p5));
	y = fq_sub(y, tmp);
	y = fq_add(y, fq_load_u(_ps_1));

	simd_fq y2 = fq_load_u(_ps_sincof_p0);
	y2 = fq_mul(y2, z);
	y2 = fq_add(y2, fq_load_u(_ps_sincof_p1));
	y2 = fq_mul(y2, z);
	y2 = fq_add(y2, fq_load_u(_ps_sincof_p2));
	y2 = fq_mul(y2, z);
	y2 = fq_mul(y2, x);
	y2 = fq_add(y2, x);

	xmm3 = poly_mask;
	y2 = fq_and(xmm3, y2);
	y = fq_and(fq_not(xmm3), y);
	y = fq_add(y, y2);

	y = fq_xor(y, sign_bit);

	return y;
}





static inline simd_fq simd_fq_cos(simd_fq x)
{
	simd_fq xmm1;
	simd_fq xmm2 = fq_zeroes();
	simd_fq xmm3;
	simd_fq y;

	simd_iq emm0;
	simd_iq emm2;

	x = fq_and(x, iq_as_fq(iq_load_u(_ps_inv_sign_mask)));

	y = fq_mul(x, fq_load_u(_ps_cephes_FOPI));

	emm2 = fq_truncate_iq(y);

	emm2 = iq_add(emm2, iq_load_u(_pi32_1));
	emm2 = iq_and(emm2, iq_load_u(_pi32_inv1));
	y = iq_convert_fq(emm2);

	emm2 = iq_sub(emm2, iq_load_u(_pi32_2));

	simd_iq tmp_invert = uq_as_iq(uq_not(iq_as_uq(emm2)));
	emm0 = iq_and(tmp_invert, iq_load_u(_pi32_4));
	emm0 = uq_as_iq(uq_lshift_29(iq_as_uq(emm0)));

	emm2 = iq_and(emm2, iq_load_u(_pi32_2));
	emm2 = iq_cmpeq(emm2, iq_zeroes());

	simd_fq sign_bit = iq_as_fq(emm0);
	simd_fq poly_mask = iq_as_fq(emm2);

	xmm1 = fq_load_u(_ps_minus_cephes_DP1);
	xmm2 = fq_load_u(_ps_minus_cephes_DP2);
	xmm3 = fq_load_u(_ps_minus_cephes_DP3);
	xmm1 = fq_mul(y, xmm1);
	xmm2 = fq_mul(y, xmm2);
	xmm3 = fq_mul(y, xmm3);
	x = fq_add(x, xmm1);
	x = fq_add(x, xmm2);
	x = fq_add(x, xmm3);

	y = fq_load_u(_ps_coscof_p0);
	simd_fq z = fq_mul(x, x);

	y = fq_mul(y, z);
	y = fq_add(y, fq_load_u(_ps_coscof_p1));
	y = fq_mul(y, z);
	y = fq_add(y, fq_load_u(_ps_coscof_p2));
	//Yes, twice
	y = fq_mul(y, z);
	y = fq_mul(y, z);
	simd_fq tmp = fq_mul(z, fq_load_u(_ps_0p5));
	y = fq_sub(y, tmp);
	y = fq_add(y, fq_load_u(_ps_1));

	simd_fq y2 = fq_load_u(_ps_sincof_p0);
	y2 = fq_mul(y2, z);
	y2 = fq_add(y2, fq_load_u(_ps_sincof_p1));
	y2 = fq_mul(y2, z);
	y2 = fq_add(y2, fq_load_u(_ps_sincof_p2));
	y2 = fq_mul(y2, z);
	y2 = fq_mul(y2, x);
	y2 = fq_add(y2, x);

	xmm3 = poly_mask;
	y2 = fq_and(xmm3, y2);
	y = fq_and(fq_not(xmm3), y);
	y = fq_add(y, y2);

	y = fq_xor(y, sign_bit);

	return y;
}






/*
sse_mathfun_extension.h - zlib license
Written by Tolga Mizrak 2016
Extension of sse_mathfun.h, which is written by Julien Pommier
Based on the corresponding algorithms of the cephes math library
This is written as an extension to sse_mathfun.h instead of modifying it, just because I didn't want
to maintain a modified version of the original library. This way switching to a newer version of the
library won't be a hassle.
Note that non SSE2 implementations of tan_ps, atan_ps, cot_ps and atan2_ps are not implemented yet.
As such, currently you need to #define USE_SSE2 to compile.
With tan_ps, cot_ps you get good precision on input ranges that are further away from the domain
borders (-PI/2, PI/2 for tan and 0, 1 for cot). See the results on the deviations for these
functions on my machine:
checking tan on [-0.25*Pi, 0.25*Pi]
max deviation from tanf(x): 1.19209e-07 at 0.250000006957*Pi, max deviation from cephes_tan(x):
5.96046e-08
   ->> precision OK for the tan_ps <<-
checking tan on [-0.49*Pi, 0.49*Pi]
max deviation from tanf(x): 3.8147e-06 at -0.490000009841*Pi, max deviation from cephes_tan(x):
9.53674e-07
   ->> precision OK for the tan_ps <<-
checking cot on [0.2*Pi, 0.7*Pi]
max deviation from cotf(x): 1.19209e-07 at 0.204303119606*Pi, max deviation from cephes_cot(x):
1.19209e-07
   ->> precision OK for the cot_ps <<-
checking cot on [0.01*Pi, 0.99*Pi]
max deviation from cotf(x): 3.8147e-06 at 0.987876517942*Pi, max deviation from cephes_cot(x):
9.53674e-07
   ->> precision OK for the cot_ps <<-
With atan_ps and atan2_ps you get pretty good precision, atan_ps max deviation is < 2e-7 and
atan2_ps max deviation is < 2.5e-7
*/

/* Copyright (C) 2016 Tolga Mizrak
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.
  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:
  1. The origin of this software must not be misrepresented; you must not
	 claim that you wrote the original software. If you use this software
	 in a product, an acknowledgment in the product documentation would be
	 appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
	 misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
  (this is the zlib license)
*/


//Mathfun Extensions

_PS_CONST( 0, 0 );
_PS_CONST( 2, 2 );
_PI32_CONST( neg1, 1 );

_PS_CONST( tancof_p0, 9.38540185543E-3f );
_PS_CONST( tancof_p1, 3.11992232697E-3f );
_PS_CONST( tancof_p2, 2.44301354525E-2f );
_PS_CONST( tancof_p3, 5.34112807005E-2f );
_PS_CONST( tancof_p4, 1.33387994085E-1f );
_PS_CONST( tancof_p5, 3.33331568548E-1f );

_PS_CONST( tancot_eps, 1.0e-4f );


static inline simd_fq simd_fq_tan(simd_fq x)
{
	simd_fq xmm1;
	simd_fq xmm2 = fq_zeroes();
	simd_fq xmm3;
	simd_fq sign_bit;
	simd_fq y;

	simd_iq emm2;

	sign_bit = x;

	x = fq_and(sign_bit, iq_as_fq(iq_load_u(_ps_inv_sign_mask)));
	sign_bit = fq_and(sign_bit, iq_as_fq(iq_load_u(_ps_sign_mask)));

	y = fq_mul(x, fq_load_u(_ps_cephes_FOPI));

	emm2 = fq_truncate_iq(y);

	emm2 = iq_add(emm2, iq_load_u(_pi32_1));
	emm2 = iq_and(emm2, iq_load_u(_pi32_inv1));
	y = iq_convert_fq(emm2);

	emm2 = iq_and(emm2, iq_load_u(_pi32_2));
	emm2 = iq_cmpeq(emm2, iq_zeroes());

	simd_fq poly_mask = iq_as_fq(emm2);

	xmm1 = fq_load_u(_ps_minus_cephes_DP1);
	xmm2 = fq_load_u(_ps_minus_cephes_DP2);
	xmm3 = fq_load_u(_ps_minus_cephes_DP3);
	xmm1 = fq_mul(y, xmm1);
	xmm2 = fq_mul(y, xmm2);
	xmm3 = fq_mul(y, xmm3);
	simd_fq z = fq_add(x, xmm1);
	z = fq_add(z, xmm2);
	z = fq_add(z, xmm3);

	simd_fq zz = fq_mul(z, z);

	y = fq_load_u(_ps_tancof_p0);
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p1));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p2));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p3));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p4));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p5));
	y = fq_mul(y, zz);
	y = fq_mul(y, z);
	y = fq_add(y, z);

		simd_fq y2 = fq_div(fq_load_u(_ps_1), y);
		y2 = fq_xor(y2, iq_as_fq(iq_load_u(_ps_sign_mask)));

	xmm3 = poly_mask;
	y = fq_and(xmm3, y);
	y2 = fq_and(fq_not(xmm3), y2);
	y = fq_or(y, y2);

	y = fq_xor(y, sign_bit);

	return y;
}

static inline simd_fq simd_fq_cot(simd_fq x)
{
	simd_fq xmm1;
	simd_fq xmm2 = fq_zeroes();
	simd_fq xmm3;
	simd_fq sign_bit;
	simd_fq y;

	simd_iq emm2;

	sign_bit = x;

	x = fq_and(sign_bit, iq_as_fq(iq_load_u(_ps_inv_sign_mask)));
	sign_bit = fq_and(sign_bit, iq_as_fq(iq_load_u(_ps_sign_mask)));

	y = fq_mul(x, fq_load_u(_ps_cephes_FOPI));

	emm2 = fq_truncate_iq(y);

	emm2 = iq_add(emm2, iq_load_u(_pi32_1));
	emm2 = iq_and(emm2, iq_load_u(_pi32_inv1));
	y = iq_convert_fq(emm2);

	emm2 = iq_and(emm2, iq_load_u(_pi32_2));
	emm2 = iq_cmpeq(emm2, iq_zeroes());

	simd_fq poly_mask = iq_as_fq(emm2);

	xmm1 = fq_load_u(_ps_minus_cephes_DP1);
	xmm2 = fq_load_u(_ps_minus_cephes_DP2);
	xmm3 = fq_load_u(_ps_minus_cephes_DP3);
	xmm1 = fq_mul(y, xmm1);
	xmm2 = fq_mul(y, xmm2);
	xmm3 = fq_mul(y, xmm3);
	simd_fq z = fq_add(x, xmm1);
	z = fq_add(z, xmm2);
	z = fq_add(z, xmm3);

	simd_fq zz = fq_mul(z, z);

	y = fq_load_u(_ps_tancof_p0);
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p1));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p2));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p3));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p4));
	y = fq_mul(y, zz);
	y = fq_add(y, fq_load_u(_ps_tancof_p5));
	y = fq_mul(y, zz);
	y = fq_mul(y, z);
	y = fq_add(y, z);

		simd_fq y2 = fq_xor(y, iq_as_fq(iq_load_u(_ps_sign_mask)));
		y = fq_div(fq_load_u(_ps_1), y);

	xmm3 = poly_mask;
	y = fq_and(xmm3, y);
	y2 = fq_and(fq_not(xmm3), y2);
	y = fq_or(y, y2);

	y = fq_xor(y, sign_bit);

	return y;
}








_PS_CONST( atanrange_hi, 2.414213562373095f );
_PS_CONST( atanrange_lo, 0.4142135623730950f );
const float PIF = 3.141592653589793238f;
const float PIO2F = 1.5707963267948966192f;
_PS_CONST( cephes_PIF, 3.141592653589793238f );
_PS_CONST( cephes_PIO2F, 1.5707963267948966192f );
_PS_CONST( cephes_PIO4F, 0.7853981633974483096f );

_PS_CONST( atancof_p0, 8.05374449538e-2f );
_PS_CONST( atancof_p1, 1.38776856032E-1f );
_PS_CONST( atancof_p2, 1.99777106478E-1f );
_PS_CONST( atancof_p3, 3.33329491539E-1f );



static inline simd_fq simd_fq_atan(simd_fq x)
{
	simd_fq sign_bit;
	simd_fq y;

	sign_bit = x;

	x = fq_and(x, iq_as_fq(iq_load_u(_ps_inv_sign_mask)));
	sign_bit = fq_and(sign_bit, iq_as_fq(iq_load_u(_ps_sign_mask)));

	simd_fq cmp0 = fq_cmpgt(x, fq_load_u(_ps_atanrange_hi));
	simd_fq cmp1 = fq_cmpgt(x, fq_load_u(_ps_atanrange_lo));
	simd_fq cmp2 = fq_and(fq_not(cmp0), cmp1);

	simd_fq y0 = fq_and(cmp0, fq_load_u(_ps_cephes_PIO2F));
	simd_fq x0 = fq_div(fq_load_u(_ps_1), x);
	x0 = fq_xor(x0, iq_as_fq(iq_load_u(_ps_sign_mask)));

	simd_fq y1 = fq_and(cmp2, fq_load_u(_ps_cephes_PIO4F));

	simd_fq x1_o = fq_sub(x, fq_load_u(_ps_1));
	simd_fq x1_u = fq_add(x, fq_load_u(_ps_1));
	simd_fq x1 = fq_div(x1_o, x1_u);

	simd_fq x2 = fq_and(cmp2, x1);
	x0 = fq_and(cmp0, x0);
	x2 = fq_or(x2, x0);
	cmp1 = fq_or(cmp0, cmp2);
	x2 = fq_and(cmp1, x2);
	x = fq_and(fq_not(cmp1), x);
	x = fq_or(x2, x);

	y = fq_or(y0, y1);


	simd_fq zz = fq_mul(x, x);
	simd_fq acc = fq_load_u(_ps_atancof_p0);
	acc = fq_mul(acc, zz);
	acc = fq_sub(acc, fq_load_u(_ps_atancof_p1));

	acc = fq_mul(acc, zz);
	acc = fq_add(acc, fq_load_u(_ps_atancof_p2));

	acc = fq_mul(acc, zz);
	acc = fq_sub(acc, fq_load_u(_ps_atancof_p3));

	acc = fq_mul(acc, zz);
	acc = fq_mul(acc, x);
	acc = fq_add(acc, x);
	y = fq_add(y, acc);

	y = fq_xor(y, sign_bit);

	return y;
}



static inline simd_fq simd_fq_atan2(simd_fq y, simd_fq x)
{
	simd_fq x_eq_0 = fq_cmpeq(x, fq_load_u(_ps_0));
	simd_fq x_gt_0 = fq_cmpgt(x, fq_load_u(_ps_0));
	simd_fq x_le_0 = fq_cmple(x, fq_load_u(_ps_0));
	simd_fq y_eq_0 = fq_cmpeq(y, fq_load_u(_ps_0));
	simd_fq x_lt_0 = fq_cmplt(x, fq_load_u(_ps_0));
	simd_fq y_lt_0 = fq_cmplt(y, fq_load_u(_ps_0));

	simd_fq zero_mask = fq_and(x_eq_0, y_eq_0);
	simd_fq zero_mask_other_case = fq_and(y_eq_0, x_gt_0);
	zero_mask = fq_or(zero_mask, zero_mask_other_case);

	simd_fq pio2_mask = fq_and(fq_not(y_eq_0), x_eq_0);
	simd_fq pio2_mask_sign = fq_and(y_lt_0, iq_as_fq(iq_load_u(_ps_sign_mask)));
	simd_fq pio2_result = fq_load_u(_ps_cephes_PIO2F);
	pio2_result = fq_xor(pio2_result, pio2_mask_sign);
	pio2_result = fq_and(pio2_mask, pio2_result);

	simd_fq pi_mask = fq_and(y_eq_0, x_le_0);
	simd_fq pi = fq_load_u(_ps_cephes_PIF);
	simd_fq pi_result = fq_and(pi_mask, pi);

	simd_fq swap_sign_mask_offset = fq_and(x_lt_0, y_lt_0);
	swap_sign_mask_offset = fq_and(swap_sign_mask_offset, iq_as_fq(iq_load_u(_ps_sign_mask)));

	simd_fq offset0 = fq_zeroes();
	simd_fq offset1 = fq_load_u(_ps_cephes_PIF);
	offset1 = fq_xor(offset1, swap_sign_mask_offset);
	
	simd_fq offset = fq_and(fq_not(x_lt_0), offset0);
	offset = fq_and(x_lt_0, offset1);

	simd_fq arg = fq_div(y, x);
	simd_fq atan_result = simd_fq_atan(arg);
	atan_result = fq_add(atan_result, offset);

	simd_fq result = fq_and(fq_not(zero_mask), pio2_result);
	//Twice in a row...is this right?
	atan_result = fq_and(fq_not(pio2_mask), atan_result);
	atan_result = fq_and(fq_not(pio2_mask), atan_result);
	result = fq_or(result, atan_result);
	result = fq_or(result, pi_result);

	return result;
}