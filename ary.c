#include "ary.h"
#include <math.h> // fabs(), isinf(), isnan()
#include <assert.h> // assert()
#include <stdio.h> // NULL

// ------------------- UTILS -------------------

// returns the minimum of a and b (or the other argument if one is NAN)
double min(double a, double b) {
	if(isnan(a)) return b;
	if(isnan(b)) return a;
	if(a < b) return a;
	return b;
}
// returns the maximum of a and b (or the other argument if one is NAN)
double max(double a, double b) {
	if(isnan(a)) return b;
	if(isnan(b)) return a;
	if(a < b) return b;
	return a;
}

// swaps values of a and b
// Requirements: a and b are not NULL
void swap(wartosc* a, wartosc* b) {
	assert(a != NULL && b != NULL);

	wartosc c = *a;
	*a = *b;
	*b = c;
}

// epsilon; the smallest positive value
const double EPS = 1e-10;
// is a equal to b (with epsilon approximation)
// or false if any of the arguments is NAN
bool eq(double a, double b) {
	return fabs(a - b) < EPS;
}
// is a less or equal to b (with epsilon approximation)
// or false if any of the arguments is NAN
bool leq(double a, double b) {
	return a < b || eq(a, b);
}
// is a greater or equal to b (with epsilon approximation)
// or false if any of the arguments is NAN
bool geq(double a, double b) {
	return a > b || eq(a, b);
}
// returns sign of a (with epsilon approximation) or 0 if a is NAN
int sgn(double a) {
	if(eq(a, 0.0) || isnan(a)) return 0;
	if(a < 0.0) return -1;
	return 1;
}
// is x equal to signed inf (positive if sign > 0, negative if sign < 0)
// Requirements: sign != 0
bool is_inf(double x, int sign) {
	assert(sign != 0);

	if(!isinf(x)) return false;
	if(x < 0.0 && sign < 0) return true;
	if(x > 0.0 && sign > 0) return true;
	return false;
}

// ------------------- CONSTRUCTORS -------------------

wartosc wartosc_dokladnosc(double x, double p) {
	assert(p > 0);

	double a = x * (100.0 - p) / 100.0;
	double b = x * (100.0 + p) / 100.0;
	return (wartosc){.first = min(a, b), .second = max(a, b), .is_flipped = false};
}

wartosc wartosc_od_do(double x, double y) {
	assert(x <= y);

	return (wartosc){.first = x, .second = y, .is_flipped = false};
}

wartosc wartosc_dokladna(double x) {
	return (wartosc){.first = x, .second = x, .is_flipped = false};
}

// ------------------- QUERIES -------------------

bool in_wartosc(wartosc w, double x) {
	if(isnan(w.first)) return false;

	if(w.is_flipped) {
		return geq(x, w.first) || leq(x, w.second);
	}
	return geq(x, w.first) && leq(x, w.second);
}

double min_wartosc(wartosc w) {
	if(isnan(w.first)) return NAN;

	// if w is flipped then it also 'contains' -inf
	if(w.is_flipped || is_inf(w.first, -1)) {
		return -HUGE_VAL;
	}
	return w.first;
}
double max_wartosc(wartosc w) {
	if(isnan(w.first)) return NAN;

	// if w is flipped then it also 'contains' +inf
	if(w.is_flipped || is_inf(w.second, 1)) {
		return HUGE_VAL;
	}
	return w.second;
}
double sr_wartosc(wartosc w) {
	double maxd = max_wartosc(w);
	double mind = min_wartosc(w);
	if(is_inf(maxd, 1) && is_inf(mind, -1)) {
		return NAN;
	}
	// note that if maxd or mind is NAN, then the return value will also be NAN - which is intended behaviour
	return (maxd + mind) / 2.0;
}

// ------------------- OPERATIONS -------------------

// returns the negation of w
// {x | -x in w}
wartosc negative(wartosc w) {
	if(isnan(w.first)) { // to preserve the invariant [*0]
		return (wartosc){.first = NAN, .second = NAN, .is_flipped = false};
	}
	return (wartosc){.first = -w.second, .second = -w.first, .is_flipped = w.is_flipped};
}

// are all endpoints of w negative
bool is_all_negative(wartosc w) {
	return leq(w.first, 0.0) && leq(w.second, 0.0);
}

wartosc plus(wartosc a, wartosc b) {
	// if both a and b are flipped, then every number can be obtained by addition
	if(a.is_flipped && b.is_flipped) {
		return (wartosc){.first = -HUGE_VAL, .second = HUGE_VAL, .is_flipped = false};
	}
	double first = a.first + b.first, second = a.second + b.second;

	// if one is flipped, then for certain arguments the result might be [-inf; +inf]
	if((a.is_flipped || b.is_flipped) && leq(first, second)) {
		return (wartosc){.first = -HUGE_VAL, .second = HUGE_VAL, .is_flipped = false};
	}
	return (wartosc){.first = first, .second = second, .is_flipped = a.is_flipped || b.is_flipped};
}
wartosc minus(wartosc a, wartosc b) {
	return plus(a, negative(b));
}

// returns the inverse of w
// {x | 1/x in w}
wartosc inverse(wartosc w) {
	// division by exactly 0.0 is undefined
	if(eq(w.first, 0.0) && eq(w.second, 0.0)) {
		return (wartosc){.first = NAN, .second = NAN, .is_flipped = false};
	}
	wartosc res = {.first = 1.0 / w.second, .second = 1.0 / w.first, .is_flipped = w.is_flipped};
	if(sgn(w.first) * sgn(w.second) == -1) {
		res.is_flipped = !w.is_flipped;
		if(eq(res.first, res.second)) { // if the segment was flipped but the endpoints are now the same
			res.first = -HUGE_VAL;
			res.second = HUGE_VAL;
			res.is_flipped = false;
		}
	}
	// handle the 'almost 0' scenarios (if one occurs, the is_flipped cancels out)
	if(eq(w.first, 0.0)) {
		res.second = HUGE_VAL;
		res.is_flipped = false;
	}
	if(eq(w.second, 0.0)) {
		res.first = -HUGE_VAL;
		res.is_flipped = false;
	}

	return res;
}

// multiply a and b if none are flipped
// Requirements: neither a nor b are flipped
wartosc mult_not_flipped(wartosc a, wartosc b) {
	assert(!a.is_flipped && !b.is_flipped);

	// mind = min{a.i * b.j | i, j in {first, second}}
	double mind	= NAN;
	mind = min(mind, a.first * b.first);
	mind = min(mind, a.first * b.second);
	mind = min(mind, a.second * b.first);
	mind = min(mind, a.second * b.second);
	// maxd = max{a.i * b.j | i, j in {first, second}}
	double maxd	= NAN;
	maxd = max(maxd, a.first * b.first);
	maxd = max(maxd, a.first * b.second);
	maxd = max(maxd, a.second * b.first);
	maxd = max(maxd, a.second * b.second);

	return (wartosc){.first = mind, .second = maxd, .is_flipped = false};
}
// multiply a and b if only one is flipped
// Requirements: *either* a or b is flipped
wartosc mult_one_flipped(wartosc a, wartosc b) {
	assert(a.is_flipped ^ b.is_flipped);

	if(a.is_flipped) swap(&a, &b); // now a is not flipped and b is 
	
	// the formulas change when working with negative-only segments,
	// so instead of writing them down, we just negate the said segments
	// and then negate the result if needed
	int sign = 1;
	if(is_all_negative(a)) {
		a = negative(a);
		sign *= -1;
	}
	if(is_all_negative(b)) {
		b = negative(b);
		sign *= -1;
	}

	double res1 = min(a.first * b.first, a.second * b.first); // first new endpoint
	double res2 = max(a.first * b.second, a.second * b.second); // second new endpoint

	if(leq(res1, res2)) { // the new segment overlaps with itself -> it is [-inf, +inf]
		return (wartosc){.first = -HUGE_VAL, .second = HUGE_VAL, .is_flipped = false};
	}

	wartosc res = (wartosc){.first = res1, .second = res2, .is_flipped = true};
	if(sign < 0) res = negative(res);
	return res;
}
// multiply a and b if both are flipped
// Requirements: both a and b are flipped
wartosc mult_both_flipped(wartosc a, wartosc b) {
	assert(a.is_flipped && b.is_flipped);

	if(in_wartosc(a, 0.0) || in_wartosc(b, 0.0)) { // if any of the segments contains 0.0, then every number can be obtained
		return (wartosc){.first = -HUGE_VAL, .second = HUGE_VAL, .is_flipped = false};
	}

	return (wartosc){
		.first = min(a.first * b.first, a.second * b.second),
		.second = max(a.first * b.second, a.second * b.first),
		.is_flipped = true
	};
}

wartosc razy(wartosc a, wartosc b) {
	if(isnan(a.first) || isnan(b.first)) { // if any of the segments is NAN, then the result is also NAN
		return (wartosc){.first = NAN, .second = NAN, .is_flipped = false};
	}
	if((eq(a.first, 0.0) && eq(a.second, 0)) || (eq(b.first, 0.0) && eq(b.second, 0))) { // special case for multiplying by [0.0, 0.0] [*1]
		return (wartosc){.first = 0.0, .second = 0.0, .is_flipped = false};
	}
	if((is_inf(a.first, -1) && is_inf(a.second, 1)) || (is_inf(b.first, -1) && is_inf(b.second, 1))) { // special case for multiplying by [-inf, inf] [*2]
		return (wartosc){.first = -HUGE_VAL, .second = HUGE_VAL, .is_flipped = false};
	}
	if(a.is_flipped && b.is_flipped) {
		return mult_both_flipped(a, b);
	}
	if(a.is_flipped || b.is_flipped) {
		return mult_one_flipped(a, b);
	}
	return mult_not_flipped(a, b);
}
wartosc podzielic(wartosc a, wartosc b) {
	return razy(a, inverse(b));
}


// ------------------- NOTES -------------------
// In the struct wartosc, the fields have the following meanings:
// - if .first = NAN, wartosc is an empty set
// - otherwise:
// - if .is_flipped = true:
//		- the set is [-inf, .second] u [.first, inf] and .second < .first holds
// - if .is_flipped = false:
//		- the set is [.first, .second] and .first <= .second holds
//
// [*0] The invariant is always preserved that wartosc A is an empty set if and only if A.first = NAN.
// Moreover, a situation where A = [not NAN, NAN] never occurs.
// Same for: [-inf, -inf], [inf, inf], [inf, -inf], [-inf, inf, is_flipped = true], [x, x, is_flipped = true]
//
// [*1] e.g. not to have NAN for test: [0, 0] * [-inf, inf] = [0, 0]
// [*2] e.g. not to have NAN for test: [-inf, inf] * [1, 0] = [-inf, inf]
//