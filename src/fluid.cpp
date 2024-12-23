#include <bits/stdc++.h>
#define vc vector

using namespace std;
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

template <int n, int k>
class Fixed {
public:
    constexpr Fixed <n, k> (int v): v(v << k) {}
    constexpr Fixed <n, k> (float f): v(f * (1 << k)) {}
    constexpr Fixed <n, k> (double f): v(f * (1 << k)) {}
    constexpr Fixed(): v(0) {}

    static constexpr Fixed from_raw(int64_t x) {
        Fixed ret;
        ret.v = x;
        return ret;
    }

    int64_t v;

    auto operator <=> (const Fixed&) const = default;
    bool operator==(const Fixed&) const = default;
};

//static constexpr Fixed inf = Fixed <32, 16>::from_raw(std::numeric_limits<int64_t>::max());
//static constexpr Fixed eps = Fixed <32, 16>::from_raw(deltas.size());

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> operator+(Fixed <n1, k1> a, Fixed <n2, k2> b) {
    if (k1 < k2) a.v += b.v >> (k2 - k1);
    else a.v += b.v << (k1 - k2);
    return a;
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> operator-(Fixed <n1, k1> a, Fixed <n2, k2> b) {
    if (k1 < k2) a.v -= b.v >> (k2 - k1);
    else a.v -= b.v << (k1 - k2);
    return a;
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> operator*(Fixed <n1, k1> a, Fixed <n2, k2> b) {
    Fixed <n1, k1> res;
    res.v = a.v;
    res.v *= b.v;
    res.v >>= k2;
    return res;
}

template <int n1, int k1>
double operator*(double a, Fixed <n1, k1> b) {
    a *= b.v>>k1;
    return a;
}

template <int n1, int k1>
Fixed <n1, k1> operator*(Fixed <n1, k1> a, double b) {
    a.v *= int(b);
    return a;
}

template <int n1, int k1>
Fixed <n1, k1> operator/(Fixed <n1, k1> a, int b) {
    return a/ Fixed <n1, k1> (b);
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> operator/(Fixed <n1, k1> a, Fixed <n2, k2> b) {
    Fixed <n1, k1> res;
    res.v = a.v;
    res.v <<= k2;
    res.v /= b.v;
    return res;
}

template <int n1, int k1>
double operator/(double a, Fixed <n1, k1> b) {
    a /= b.v>>k1;
    return a;
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> &operator+=(Fixed <n1, k1> &a, Fixed <n2, k2> b) {
    return a = a + b;
}

template <int n1, int k1>
Fixed <n1, k1> &operator+=(Fixed <n1, k1> &a, double b) {
    a.v += int(b) << k1;
    return a;
}

template <int n1, int k1>
double &operator+=(double& a, Fixed <n1, k1> &b) {
    a += b.v>>k1;
    return a;
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> &operator-=(Fixed <n1, k1> &a, Fixed <n2, k2> b) {
    return a = a - b;
}

template <int n1, int k1>
Fixed <n1, k1> &operator-=(Fixed <n1, k1> &a, double b) {
    Fixed <n1, k1> b_new = b;
    return a = a - b_new;
}

template <int n1, int k1>
double &operator-=(double &a, Fixed <n1, k1> &b) {
    a -= b.v>>k1;
    return a;
}

template <int n, int k>
Fixed <n, k> &operator*=(Fixed <n, k> &a, double b) {
    return a = a * Fixed <n, k> (b);
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> &operator*=(Fixed <n1, k1> &a, Fixed <n2, k2> b) {
    return a = a * b;
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> &operator/=(Fixed <n1, k1> &a, Fixed <n2, k2> b) {
    return a = a / b;
}

template <int n, int k>
Fixed <n, k> operator-(Fixed <n, k> x) {
    return Fixed <n, k>::from_raw(-x.v);
}

template <int n, int k>
Fixed <n, k> abs(Fixed <n, k> x) {
    if (x.v < 0) {
        x.v = -x.v;
    }
    return x;
}

template <int n, int k>
ostream &operator<<(ostream &out, Fixed <n, k> x) {
    return out << x.v / (double) (1 << k);
}

int UT = 0;

mt19937 rnd(1337);

template <int n, int k>
Fixed <n, k> random01() {
    return Fixed <n, k>::from_raw((rnd() & ((1 << k) - 1)));
}

//start
ifstream fin("input.txt");
template <int n, int k>
Fixed <n, k> read_new() {
    float a;
    fin >> a;
    return Fixed <n, k> (a);
}

template <typename p_type, typename v_type, typename v_flow_type, size_t N, size_t M>
class simulation {
public:
    int dirs[N][M]{};
    p_type rho[256];
    p_type g = 0.5;
    p_type p[N][M]{};
    p_type old_p[N][M];
    int last_use[N][M]{};

    template<typename T>
    struct VectorField {
        array<T, deltas.size()> v[size_t(N)][M];
        T &add(int x, int y, int dx, int dy, T dv) {
            return get(x, y, dx, dy) += dv;
        }

        T &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    VectorField <v_type> velocity{};
    VectorField <v_flow_type> velocity_flow{};

    char field[N][M + 1];

    tuple<p_type, bool, pair<int, int>> propagate_flow(int x, int y, p_type lim) {
        last_use[x][y] = UT - 1;
        p_type ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                v_type cap = velocity.get(x, y, dx, dy);
                v_type flow = 0;
                flow += velocity_flow.get(x, y, dx, dy);

                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                v_type nlim = 0;
                nlim += lim;
                auto vp = min(nlim, cap - flow);
                if (last_use[nx][ny] == UT - 1) {
                    v_flow_type nvp = 0;
                    nvp += vp;
                    velocity_flow.add(x, y, dx, dy, nvp);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    p_type to_ret = 0;
                    to_ret += nvp;
                    return {to_ret, 1, {nx, ny}};
                }
                p_type to_ret = 0;
                to_ret += vp;
                auto [t, prop, end] = propagate_flow(nx, ny, to_ret);
                ret += t;
                if (prop) {
                    v_flow_type nt = 0;
                    nt += t;
                    velocity_flow.add(x, y, dx, dy, nt);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
                    return {t, prop && end != pair(x, y), end};
                }
            }
        }
        last_use[x][y] = UT;
        return {ret, 0, {0, 0}};
    }

    struct ParticleParams {
        char type;
        p_type cur_p;
        array<v_type, deltas.size()> v;

        void swap_with(int x, int y, simulation * ths) {
            swap(ths->field[x][y], type);
            swap(ths->p[x][y], cur_p);
            swap(ths->velocity.v[x][y], v);
        }
    };

    void propagate_stop(int x, int y, bool force = false) {
        if (!force) {
            bool stop = true;
            for (auto [dx, dy] : deltas) {
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > 0) {
                    stop = false;
                    break;
                }
            }
            if (!stop) {
                return;
            }
        }
        last_use[x][y] = UT;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > 0) {
                continue;
            }
            propagate_stop(nx, ny);
        }
    }

    p_type move_prob(int x, int y) {
        p_type sum = 0;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                continue;
            }
            auto v = velocity.get(x, y, dx, dy);
            if (v < 0) {
                continue;
            }
            sum += v;
        }
        return sum;
    }

    bool propagate_move(int x, int y, bool is_first) {
        last_use[x][y] = UT - is_first;
        bool ret = false;
        int nx = -1, ny = -1;
        do {
            array<p_type, deltas.size()> tres;
            p_type sum = 0;
            for (size_t i = 0; i < deltas.size(); ++i) {
                auto [dx, dy] = deltas[i];
                int nx = x + dx, ny = y + dy;
                if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
                    tres[i] = sum;
                    continue;
                }
                auto v = velocity.get(x, y, dx, dy);
                if (v < 0) {
                    tres[i] = sum;
                    continue;
                }
                sum += v;
                tres[i] = sum;
            }

            if (sum == 0) {
                break;
            }

            p_type p1 = 0;
            Fixed <32, 16> xx = random01 <32, 16>() * sum;
            p1 += xx;
            size_t d = std::ranges::upper_bound(tres, p1) - tres.begin();

            auto [dx, dy] = deltas[d];
            nx = x + dx;
            ny = y + dy;
            assert(velocity.get(x, y, dx, dy) > 0 && field[nx][ny] != '#' && last_use[nx][ny] < UT);

            ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
        } while (!ret);
        last_use[x][y] = UT;
        for (size_t i = 0; i < deltas.size(); ++i) {
            auto [dx, dy] = deltas[i];
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < 0) {
                propagate_stop(nx, ny);
            }
        }
        if (ret) {
            if (!is_first) {
                ParticleParams pp{};
                pp.swap_with(x, y, this);
                pp.swap_with(nx, ny, this);
                pp.swap_with(x, y, this);
            }
        }
        return ret;
    }

    void init() {
        double cur;
        fin >> cur;
        rho[' '] = cur;

        fin >> cur;
        rho['.'] = cur;

        fin >> cur;
        g = cur;

        string str;
        getline(fin, str);;
        for (int i = 0; i < N; i++) {

            getline(fin, str);
            for (int j = 0; j < M; j++) {
                //field[i][j] = '#';
                field[i][j] = str[j];
                //cout << field[i][j];
            }
            //cout << field[i] << "\n";
            //cout << "\n";
        }
    }

    void work() {
        for (size_t x = 0; x < N; ++x) {
            for (size_t y = 0; y < M; ++y) {
                if (field[x][y] == '#')
                    continue;
                for (auto [dx, dy]: deltas) {
                    dirs[x][y] += (field[x + dx][y + dy] != '#');
                }
            }
        }

        for (size_t i = 0; i < T; ++i) {
            string filename = to_string(i);
            filename += "output.txt"; //result will be saved
            //cout << filename << ' ';
            ofstream fout (filename);

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    fout << field[x][y];
                }
                fout << "\n";
            }

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    fout << p[x][y] << ' ';
                }
                fout << "\n";
            }

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    fout << last_use[x][y] << ' ';
                }
                fout << "\n";
            }

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    fout << dirs[x][y] << ' ';
                }
                fout << "\n";
            }

            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    for (size_t ii = 0; ii < deltas.size(); ++ii) {
                        auto [dx, dy] = deltas[ii];
                        int nx = x + dx, ny = y + dy;
                        auto v = velocity.get(x, y, dx, dy);
                        fout << v;
                    }
                }
            }

            p_type total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#') {
                        v_type zero = 0;
                        zero += g;
                        velocity.add(x, y, 1, 0, zero);
                    }
                }
            }

            // Apply forces from p
            memcpy(old_p, p, sizeof(p));
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy]: deltas) {
                        int nx = x + dx, ny = y + dy;
                        if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
                            auto delta_p = old_p[x][y] - old_p[nx][ny];
                            v_type force = 0;
                            force += delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);

                            v_type plus = force / dirs[x][y];
                            p[x][y] -= plus;
                            plus = force / dirs[x][y];
                            total_delta_p -= plus;
                        }
                    }
                }
            }

            // Make flow from velocities
            velocity_flow = {};
            bool prop = false;
            do {
                UT += 2;
                prop = 0;
                for (size_t x = 0; x < N; ++x) {
                    for (size_t y = 0; y < M; ++y) {
                        if (field[x][y] != '#' && last_use[x][y] != UT) {
                            auto [t, local_prop, _] = propagate_flow(x, y, 1);
                            if (t > 0) {
                                prop = 1;
                            }
                        }
                    }
                }
            } while (prop);

            // Recalculate p with kinetic energy
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    for (auto [dx, dy]: deltas) {
                        auto old_v = velocity.get(x, y, dx, dy);
                        v_type new_v = 0;
                        new_v += velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            v_type force = (old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                v_type plus = force / dirs[x][y];
                                p[x][y] += plus;
                                plus = force / dirs[x][y];
                                total_delta_p += plus;
                            } else {
                                v_type plus = force / dirs[x + dx][y + dy];
                                p[x + dx][y + dy] += plus;
                                plus = force / dirs[x + dx][y + dy];
                                total_delta_p += plus;
                            }
                        }
                    }
                }
            }

            UT += 2;
            prop = false;
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] != '#' && last_use[x][y] != UT) {
                        p_type ff = 0;
                        Fixed <32, 16> xx = random01<32, 16>();
                        ff += xx;
                        if (ff < move_prob(x, y)) {
                            prop = true;
                            propagate_move(x, y, true);
                        } else {
                            propagate_stop(x, y, true);
                        }
                    }
                }
            }

            if (prop) {
                cout << "Tick " << i << ":\n";
                for (size_t x = 0; x < N; ++x) {
                    //cout << field[x] << "\n";
                }
            }
        }
    }
};

int main() {
    simulation <Fixed <32, 16>, Fixed <32, 16>, Fixed <32, 16>, 36, 84> x; //you can choose types and sizes
    x.init();
    x.work();
    return 0;
}
