#include <bits/stdc++.h>
#define vc vector

using namespace std;
// constexpr size_t N = 14, M = 5;
constexpr size_t T = 1'000'000;
constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

// char field[N][M + 1] = {
//     "#####",
//     "#.  #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#.# #",
//     "#...#",
//     "#####",
//     "#   #",
//     "#   #",
//     "#   #",
//     "#####",
// };

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

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> &operator+=(Fixed <n1, k1> &a, Fixed <n2, k2> b) {
    return a = a + b;
}

template <int n1, int n2, int k1, int k2>
Fixed <n1, k1> &operator-=(Fixed <n1, k1> &a, Fixed <n2, k2> b) {
    return a = a - b;
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

template <int n1, int n2, int k1, int k2>
Fixed <n2, k2> to(Fixed <n1, k1> x) {
    Fixed <n2, k2> res;
    res.v = x.v;
    if (k1 < k2) res.v <<= (k2 - k1);
    else res.v >>= (k1 - k2);
    return res;
}



int UT = 0;


mt19937 rnd(1337);

template <int n, int k>
Fixed <n, k> random01() {
    return Fixed <n, k>::from_raw((rnd() & ((1 << k) - 1)));
}

//start
ifstream fin("D:\\homework\\input.txt");
template <int n, int k>
Fixed <n, k> read_new() {
    float a;
    fin >> a;
    return Fixed <n, k> (a);
}
//#define T_r Fixed <30, 13>
//#define T_def  Fixed <32, 16>

constexpr size_t N = 36, M = 84;
template <typename T_r, typename T_def>
class simulation {
public:
    int dirs[N][M]{};
    T_r rho[256];
    T_r g = 0.5;
    T_r p[N][M]{};
    T_r old_p[N][M];
    int last_use[N][M]{};

    struct VectorField {
        array<Fixed <30, 13>, deltas.size()> v[size_t(N)][M];
        Fixed <30, 13> &add(int x, int y, int dx, int dy, Fixed <30, 13> dv) {
            return get(x, y, dx, dy) += dv;
        }

        Fixed <30, 13> &get(int x, int y, int dx, int dy) {
            size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
            assert(i < deltas.size());
            return v[x][y][i];
        }
    };

    VectorField velocity{}, velocity_flow{};

    char field[N][M + 1];

    tuple<T_r, bool, pair<int, int>> propagate_flow(int x, int y, T_r lim) {
        last_use[x][y] = UT - 1;
        T_r ret = 0;
        for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
                auto cap = velocity.get(x, y, dx, dy);
                auto flow = velocity_flow.get(x, y, dx, dy);
                if (flow == cap) {
                    continue;
                }
                // assert(v >= velocity_flow.get(x, y, dx, dy));
                auto vp = min(lim, cap - flow);
                if (last_use[nx][ny] == UT - 1) {
                    velocity_flow.add(x, y, dx, dy, vp);
                    last_use[x][y] = UT;
                    // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
                    return {vp, 1, {nx, ny}};
                }
                auto [t, prop, end] = propagate_flow(nx, ny, vp);
                ret += t;
                if (prop) {
                    velocity_flow.add(x, y, dx, dy, t);
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
        T_r cur_p;
        array<T_r, deltas.size()> v;

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

    T_def move_prob(int x, int y) {
        T_def sum = 0;
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
            std::array<Fixed <32, 16>, deltas.size()> tres;
            Fixed <32, 16> sum = 0;
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

            Fixed p1 = random01 <32, 16>() * sum;
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
            filename += "output.txt";
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

            int v1, v2, v3, v4;
            auto v = velocity.get(v1, v2, v3, v4);
            fout << v << "\n";

            T_def total_delta_p = 0;
            // Apply external forces
            for (size_t x = 0; x < N; ++x) {
                for (size_t y = 0; y < M; ++y) {
                    if (field[x][y] == '#')
                        continue;
                    if (field[x + 1][y] != '#')
                        velocity.add(x, y, 1, 0, g);
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
                            auto force = delta_p;
                            auto &contr = velocity.get(nx, ny, -dx, -dy);
                            if (contr * rho[(int) field[nx][ny]] >= force) {
                                contr -= force / rho[(int) field[nx][ny]];
                                continue;
                            }
                            force -= contr * rho[(int) field[nx][ny]];
                            contr = 0;
                            velocity.add(x, y, dx, dy, force / rho[(int) field[x][y]]);
                            p[x][y] -= force / dirs[x][y];
                            total_delta_p -= force / dirs[x][y];
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
                        auto new_v = velocity_flow.get(x, y, dx, dy);
                        if (old_v > 0) {
                            assert(new_v <= old_v);
                            velocity.get(x, y, dx, dy) = new_v;
                            auto force = (old_v - new_v) * rho[(int) field[x][y]];
                            if (field[x][y] == '.')
                                force *= 0.8;
                            if (field[x + dx][y + dy] == '#') {
                                p[x][y] += force / dirs[x][y];
                                total_delta_p += force / dirs[x][y];
                            } else {
                                p[x + dx][y + dy] += force / dirs[x + dx][y + dy];
                                total_delta_p += force / dirs[x + dx][y + dy];
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
                        if (random01<32, 16>() < move_prob(x, y)) {
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
    simulation <Fixed <30, 13>, Fixed <32, 16>> x;
    x.init();
    x.work();
    return 0;
}
