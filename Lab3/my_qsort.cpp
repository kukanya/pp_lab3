#include "stdafx.h"

void my_qsort(int* a, int left, int right) {
	int l = left, r = right, swp;
	int pivot = a[(l+r)/2];

	while (l <= r) {
		while (a[l] < pivot) ++l;
		while (a[r] > pivot) --r;
		if (l <= r) {
			swp = a[l];	a[l] = a[r]; a[r] = swp;
			++l; --r;
		}
	}
	if (r > left) my_qsort(a, left, r);
	if (l < right) my_qsort(a, l, right);
}

void my_qsort_bounded(int* a, double* b, int left, int right) {
	int l = left, r = right, swp;
	double dswp;
	int pivot = a[(l+r)/2];

	while (l <= r) {
		while (a[l] < pivot) ++l;
		while (a[r] > pivot) --r;
		if (l <= r) {
			swp = a[l];	a[l] = a[r]; a[r] = swp;
			dswp = b[l]; b[l] = b[r]; b[r] = dswp;
			++l; --r;
		}
	}
	if (r > left) my_qsort_bounded(a, b, left, r);
	if (l < right) my_qsort_bounded(a, b, l, right);
}