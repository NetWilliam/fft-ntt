#include <algorithm>
#include <complex>
#include <iostream>
#include <vector>

using namespace std;

const long double PI = 3.141592653589793238L;

using complex_type = complex<long double>;

int size_enlarge2power_of2(size_t len)
{
    if ((len ^ (len - 1)) == 0) {
        return len;
    }
    return len + (~len) + 1;
}

template <typename T>
void butterfly_exchange(T* array, size_t len)
{
    //size_t len = array.size();
    for (int i = 1, j = len / 2; i < len - 1; ++i) {
        if (i < j)
            swap(array[i], array[j]);
        int k = len / 2;
        while (j >= k) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

vector<int> multi(vector<int>& f, vector<int>& g)
{
    vector<int> ret;
    return ret;
}

vector<int> fft_multi(vector<int>& f, vector<int>& g)
{
    int f_enlarge_len = size_enlarge2power_of2(f.size());
    int g_enlarge_len = size_enlarge2power_of2(g.size());
    int max_len = f_enlarge_len > g_enlarge_len ? f_enlarge_len : g_enlarge_len;
    f.resize(max_len);
    g.resize(max_len);
    return multi(f, g);
}

int main()
{
    vector<int> f = { 5, 3, 9, 6 }; // f(x) = 5 + 3x + 9x^2 + 6x^3
    vector<int> g = { 6, 7, 2, 4 }; // g(x) = 6 + 7x + 2x^2 + 4x^3

    //0, 1, 2, 3, 4, 5, 6, 7 => butterfly_exchange => 0, 4, 2, 6, 1, 5, 3, 7
    vector<int> butterfly_exchange_arr(8);
    generate(butterfly_exchange_arr.begin(), butterfly_exchange_arr.end(), []() { static int i = 0; return i++; });
    vector<int> butterfly_exchange_test_arr = { 0, 4, 2, 6, 1, 5, 3, 7 };

    butterfly_exchange(butterfly_exchange_arr.data(), butterfly_exchange_arr.size());

    if (butterfly_exchange_test_arr == butterfly_exchange_arr) {
        cout << "butterfly_exchange function test passed!" << endl;
    } else {
        cout << "butterfly_exchange function test failed!" << endl;
        cout << "expect:" << endl;
        for_each(butterfly_exchange_test_arr.begin(), butterfly_exchange_test_arr.end(), [](int i) { cout << i << ','; });
        cout << "\nactually get:" << endl;
        for_each(butterfly_exchange_arr.begin(), butterfly_exchange_arr.end(), [](int i) { cout << i << ','; });
        cout << endl;
    }

    // expect: h(x) = 30 + 53x + 85x^2 + 125x^3 + 72x^4 + 48x^5 + 24x^6
    vector<int> h = fft_multi(f, g);
    vector<int> h_ground_true = { 30, 53, 85, 125, 72, 48, 24 };

    if (h == h_ground_true) {
        cout << "answer is right!" << endl;
    } else {
        cout << "answer is wrong!" << endl;
    }
}
