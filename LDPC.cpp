//#include <bits/stdc++.h>
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>
#include <map>
#include <unordered_map>
#include <mutex>
#include <limits.h>

/// Новая фича: использование unordered_map и разложение по последнему столбцу

using namespace std;

bool debug = false;
bool cheating = true;

vector<vector<int>> matrix;
vector<vector<int>> ss;
vector<thread> thread_pool;
unsigned numt = thread::hardware_concurrency();
double prop_non_zero = 0.8;
vector<vector<unsigned long long>> permsnsum;

struct VectorHasher {
    size_t operator()(const vector<int> &V) const {
        size_t hash = 0;
        for(auto &i : V) {
            hash |= (1 << i);
        }
        return hash;
    }
};

unordered_map<vector<int>, unsigned long long, VectorHasher> hashed_matrix;
mutex hmm;



unsigned long long get_matrix_perm(vector<int> &subset) {
    unsigned long long answer = ULLONG_MAX;

    hmm.lock();
    auto it = hashed_matrix.find(subset);
    if(it != hashed_matrix.end()) {
        /// Если найден ответ в памяти, то возвращаем ответ
        answer = it->second;
    }
    hmm.unlock();
    return answer;
}

void save_matrix_perm(vector<int> &subset, unsigned long long perm) {
    hmm.lock();
    hashed_matrix[subset] = perm;
    hmm.unlock();
}

void search(int k, int n, long unsigned int maxn, vector<int> &subset, vector<vector<int>> &ss) {

    if(subset.size() == maxn) {
        ss.push_back(subset);
        return;
    }
    if (k == n+1) {
        // обрабатываем подмножество
    } else {
        subset.push_back(k);
        search(k+1, n, maxn, subset, ss);
        subset.pop_back();
        search(k+1, n, maxn, subset, ss);
    }
}

unsigned long long calc_subperm(vector<int> x_arr,
                                vector<int> y_arr,
                                int del_x,
                                int del_y,
                                bool isDecompose)
{
    long unsigned int i;
    // Удаляем индексы ненужных столбцов и строк
    if( del_x != -1 || del_y != -1) {
        x_arr.erase(x_arr.begin() + del_x);
        y_arr.erase(y_arr.begin() + del_y);
    }
    if(x_arr.size() != y_arr.size()) {
        cout << "Внимание! Ты плохой программист!";
        return 0;
    }

    // Теперь считаем перманент
    unsigned long long summ=0;
    if((x_arr.size() !=2) && isDecompose) {
        for(i=0; i<x_arr.size(); ++i) {
            if(matrix[x_arr[i]][y_arr[0 + y_arr.size() - 1]] != 0)
                summ += matrix[x_arr[i]][y_arr[0 + y_arr.size() - 1]] * calc_subperm(x_arr, y_arr, i, 0 + y_arr.size() - 1, isDecompose);
        }
    } else {
        /// Считаем по-старому
        vector<int> permutation;
        unsigned long long product;

        for(i=0; i<x_arr.size(); ++i)
            permutation.push_back(i);

        do {
            product = 1;
            for(i = 0; i < permutation.size(); ++i) {
                if(matrix[x_arr[i]][y_arr[permutation[i]]] == 1) continue;
                product *= matrix[x_arr[i]][y_arr[permutation[i]]];
                if(product == 0)
                    break;
            }
            if(debug) {
                cout << "product = " << product << '\n';
            }
            summ += product;
        } while(next_permutation(permutation.begin(), permutation.end()));

    }

    return summ;
}

void perm(unsigned long index_to_delete, unsigned long ssnum) {

    vector<int> subset = ss[ssnum];
    // Удаляем один из элементов подмножества
    subset.erase(subset.begin() + index_to_delete);

    unsigned long int i, j, k;
    vector<int> x_arr(subset.size());
    for(i=0; i<subset.size(); ++i) {
        subset[i] -= 1;
        x_arr[i] = i;
    }

    /// Запрашиваем хэшированное значение и если оно есть, возвращаемся с ответом.
    permsnsum[ssnum][index_to_delete] = get_matrix_perm(subset);
    if(permsnsum[ssnum][index_to_delete] != ULLONG_MAX)
        return;


    if(debug) {
        cout << "Подмножество c удаленным столцом:\n";
        for(k=0; k<subset.size(); ++k)
            cout << subset[k] << '\t';
        cout << '\n';
    }
    permsnsum[ssnum][index_to_delete] = 0;

    /// Проверяем матрицу на единичность
    bool is_1 = true;
    if(cheating) {
        for(i=0; i<matrix.size() && is_1; ++i) {
            for(j=0; j<matrix.size() && is_1; ++j) {
                if(matrix[i][subset[j]] != 1)
                    is_1 = false;
            }
        }
    }

    /// Считаем результат по формуле факториала
    if(is_1 && cheating) {
        permsnsum[ssnum][index_to_delete] = 1;
        for(i=2; i<=matrix.size(); ++i) {
            permsnsum[ssnum][index_to_delete] *= i;
        }
        return;
    }


        permsnsum[ssnum][index_to_delete] = calc_subperm(x_arr, subset, -1, -1, true);
    //Сохраняем результаты вычислений
    save_matrix_perm(subset, permsnsum[ssnum][index_to_delete]);
}

void make_subsets(vector<vector<int>> &ss,
                  int n ,
                  int maxn) {
    vector<int> subset;

    search(1, n, maxn, subset, ss);
}

void perm_worker(int id, int num_workers) {
    unsigned long i, j;
    unsigned long long sum=0;
    for(i=id; i < ss.size(); i+= num_workers) {
        sum=0;
        for(j = 0; j < ss[i].size(); ++j) {
            perm(j, i);
            sum+= permsnsum[i][j];
        }
        permsnsum[i][ss[i].size()] = sum;

        if(i%500 == 0)
            cout << "Подмножество " << i << "/" << ss.size() << " готово \n";

    }
}

void run_perms_and_sum() {
    unsigned int i;
    unsigned long long sum = 0;
    thread_pool.clear();
    /// Запускаем потоки
    for(i=0; i<numt; ++i){
        thread_pool.push_back(thread(perm_worker, i, numt));
    }

    /// Ждем выполнения потоков
    for(i=0; i<numt; ++i) {
        thread_pool[i].join();
    }
}


int main() {
//    ios::sync_with_stdio(0);
//    cin.tie(0);

//    freopen("/home/nikolay/prog/t10.txt", "r", stdin);
    cout << "Число потоков: " << numt << "\n";

    unsigned long int n, m, i, j, k;
    cin >> m >> n;

    // Считываем матрицу массив
    matrix.resize(m);
    for ( i=0; i<m; ++i) {
        matrix[i].resize(n);
        for( j=0; j<n; ++j) {
            cin >> matrix[i][j];
        }
    }

    if(debug) {
        for ( i=0; i<matrix.size(); ++i) {
            for( j=0; j<matrix[i].size(); ++j) {
                cout << matrix[i][j] << '\t';
            }
            cout << '\n';
        }
    }

    // Создаем массив перестановок
    make_subsets(ss, n, m+1);

    if(debug) {
        cout << "Результат: \n";
        for (i=0; i<ss.size(); ++i) {
            for(j=0; j<ss[i].size(); ++j) {
                cout << ss[i][j] << '\t';
            }
            cout << '\n';
        }
    }

    unsigned long long sum = 0;
    unsigned long long minsum  = 0;

    permsnsum.resize(ss.size());
    for(i=0; i<permsnsum.size(); ++i) {
        permsnsum[i].resize(m+2);
    }

    run_perms_and_sum();

    for (i=0; i<ss.size(); ++i) {
        sum = permsnsum[i][ss[i].size()];
        if(debug) {
            cout << "sum=  " << sum << '\n';
        }

        if((minsum == 0) && (sum != 0)) {
            minsum = sum;
        }
        if(sum!=0) {
            minsum = min(sum, minsum);
        }
    }

    cout << "minsum = " << minsum << '\n';

    return 0;
}
