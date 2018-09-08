#include <algorithm>
#include <assert.h>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/range.hpp>
#include <complex>
#include <iostream>
#include <vector>

using namespace std;

// const long double PI = 3.141592653589793238L;
const long double PI = acos(-1.0);

using complex_type = complex<long double>;

int size_enlarge2power_of2(size_t len)
{
    if ((len & (len - 1)) == 0)
    {
        return len;
    }
    return len + (~len) + 1;
}

template <typename T> void butterfly_exchange(T *array, size_t len)
{
    // size_t len = array.size();
    assert((len & (len - 1)) == 0);
    assert(len != 0);
    for (int i = 1, j = len / 2; i < len - 1; ++i)
    {
        if (i < j)
            swap(array[i], array[j]);
        int k = len / 2;
        while (j >= k)
        {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

template <typename T> void fft(T *y, int len, int sign)
{
    /*
    for (int h = 2; h <= len; h <<= 1)
    {
        complex_type wn(cos(sign * 2 * PI / h), sin(sign * 2 * PI / h));
        for (int j = 0; j < len; j += h)
        {
            complex_type w(1, 0);
            for (int k = j; k < j + h / 2; ++k)
            {
                auto u = y[k];
                auto t = w * y[k + h / 2];
                y[k] = u + t;
                y[k + h / 2] = u - t;
                w = w * wn;
            }
        }
    }
    //*/
    for (int stride = 2; stride <= len; stride <<= 1)
    {
        complex_type wn(cos(sign * 2 * PI / stride), sin(sign * 2 * PI / stride));
        for (int step = 0; step < len; step += stride)
        {
            complex_type w(1, 0);
            for (int k = step; k < step + stride / 2; k += 1)
            {
                auto u = y[k];
                auto v = w * y[k + stride / 2];
                y[k] = u + v;
                y[k + stride / 2] = u - v;
                w = w * wn;
            }
        }
    }
    //*/
    if (sign == -1)
        for_each(y, y + len, [len](T &y_ele) { y_ele /= len; });
}

template <typename T> void print_c(T &c)
{
    for_each(c.begin(), c.end(), [](const typename T::value_type &v) { cout << v << ", "; });
    cout << endl;
}

vector<int> fft_multi(vector<int> &f, vector<int> &g)
{
    int f_enlarge_len = size_enlarge2power_of2(f.size());
    int g_enlarge_len = size_enlarge2power_of2(g.size());
    int max_len = f_enlarge_len > g_enlarge_len ? f_enlarge_len : g_enlarge_len;
    cout << "f.size(): " << f.size() << " g.size(): " << g.size() << " max_len: " << max_len << endl;
    f.resize(max_len);
    g.resize(max_len);

    vector<complex_type> y_f(2 * max_len);
    vector<complex_type> y_g(2 * max_len);
    transform(f.begin(), f.end(), y_f.begin(), [](int x) -> complex_type { return complex_type(x, 0); });
    transform(g.begin(), g.end(), y_g.begin(), [](int x) -> complex_type { return complex_type(x, 0); });
    print_c(y_f);
    print_c(y_g);

    cout << "fft transform, butterfly_exchanging ......, y_f.size(): " << y_f.size() << endl;
    butterfly_exchange(y_f.data(), y_f.size());
    butterfly_exchange(y_g.data(), y_g.size());
    print_c(y_f);
    print_c(y_g);

    cout << "fft transform, sampling ......" << endl;
    fft(y_f.data(), y_f.size(), 1);
    fft(y_g.data(), y_g.size(), 1);
    print_c(y_f);
    print_c(y_g);

    cout << "samples multiplying ......" << endl;
    for (int i = 0; i < y_f.size(); ++i)
    {
        y_f[i] *= y_g[i];
    }
    print_c(y_f);

    cout << "dft transform, butterfly_exchanging ......" << endl;
    butterfly_exchange(y_f.data(), y_f.size());
    cout << "dft transform, sampling ......" << endl;
    fft(y_f.data(), y_f.size(), -1);

    cout << "multi done, returning ......" << endl;
    vector<int> ret(y_f.size());
    transform(y_f.begin(), y_f.end(), ret.begin(),
              [](complex_type &c) -> int { return static_cast<int>(c.real() + 0.5); });

    print_c(y_f);
    return ret;
}

int main()
{
    vector<int> f = {5, 3, 9, 6}; // f(x) = 5 + 3x + 9x^2 + 6x^3
    vector<int> g = {6, 7, 2, 4}; // g(x) = 6 + 7x + 2x^2 + 4x^3

    // 0, 1, 2, 3, 4, 5, 6, 7 => butterfly_exchange => 0, 4, 2, 6, 1, 5, 3, 7
    vector<int> butterfly_exchange_arr(8);
    generate(butterfly_exchange_arr.begin(), butterfly_exchange_arr.end(), []() {
        static int i = 0;
        return i++;
    });
    vector<int> butterfly_exchange_test_arr = {0, 4, 2, 6, 1, 5, 3, 7};

    butterfly_exchange(butterfly_exchange_arr.data(), butterfly_exchange_arr.size());

    if (butterfly_exchange_test_arr == butterfly_exchange_arr)
    {
        cout << "butterfly_exchange function test passed!" << endl;
    }
    else
    {
        cout << "butterfly_exchange function test failed!" << endl;
        cout << "expect:" << endl;
        for_each(butterfly_exchange_test_arr.begin(), butterfly_exchange_test_arr.end(),
                 [](int i) { cout << i << ','; });
        cout << "\nactually get:" << endl;
        for_each(butterfly_exchange_arr.begin(), butterfly_exchange_arr.end(), [](int i) { cout << i << ','; });
        cout << endl;
    }

    // expect: h(x) = 30 + 53x + 85x^2 + 125x^3 + 72x^4 + 48x^5 + 24x^6
    vector<int> h = fft_multi(f, g);
    vector<int> h_ground_true = {30, 53, 85, 125, 72, 48, 24, 0};

    if (h == h_ground_true)
    {
        cout << "answer is right!" << endl;
    }
    else
    {
        cout << "answer is wrong!" << endl;
        cout << "butterfly_exchange function test failed!" << endl;
        cout << "expect:" << endl;
        for_each(h_ground_true.begin(), h_ground_true.end(), [](int i) { cout << i << ','; });
        cout << "\nactually get:" << endl;
        for_each(h.begin(), h.end(), [](int i) { cout << i << ','; });
        cout << endl;
    }
}
