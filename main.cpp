#include <iostream>
#include <fstream>
#include <chrono>

using namespace std;

void heapSort(int *&priority, int *&perm, int &n);

long algorytmNEH(int &n, int &m, int **&p, long &cmax, int *&perm);


void
modyfikacjaIR4(int &m, int **&p, int **&r, int **&q, int &zad, int &t, int &zad_kryt, long &cmax, int *&perm,
               int *&dist, int &temp);
void
modyfikacjaIRR4(int &m, int **&p, int **&r, int **&q, int &zad, int &t, int &zad_kryt, long &cmax, int *&perm,
                int *&dist, int &temp);

int main() {
    int n, m;
    long cmax;
    int **p;
    int *perm;
    int reps=1;
    chrono::duration<double> czas_obliczen_NEH = chrono::duration<double>::zero();
    chrono::duration<double> elapsed_seconds;

    string plik_wyj = "permutacje.txt";
    string plik_wej = "neh.data.txt";
    ifstream plik;
    ofstream wyjscie;
    string linia;
    wyjscie.open(plik_wyj, ofstream::trunc);
    wyjscie.close();
    plik.open(plik_wej);

    while (!plik.eof()) {
        getline(plik, linia);
        if (linia[0] == 'd' && linia[1] == 'a') {
            plik >> n >> m;

            p = new int *[m + 1];
            for (int i = 1; i <= m; ++i)
                p[i] = new int[n + 1];

            for (int j = 1; j <= n; ++j)
                for (int i = 1; i <= m; ++i)
                    plik >> p[i][j];
            perm = new int[n + 1];
            chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
            for (int i=0; i<reps; ++i)
                algorytmNEH(n, m, p, cmax, perm);
            chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
            elapsed_seconds = end - start;
            czas_obliczen_NEH += elapsed_seconds;

            wyjscie.open(plik_wyj, ofstream::app);
            // WYPISZ CO POTRZEBUJESZ
            //for (int j = 1; j <= n; ++j) wyjscie << perm[j] << " ";
            wyjscie << endl << cmax ;//<< endl;
            //wyjscie << n << " " << m << " ";
            //wyjscie << /*"Czas obliczen: " <<*/ elapsed_seconds.count() << endl;
            wyjscie.close();

            for (int i = 1; i <= m; ++i)
                delete[] p[i];
            delete[] p;
            delete[] perm;
        }
    }
    plik.close();
    wyjscie.open(plik_wyj, ofstream::app);
    //wyjscie << "Calkowity czas obliczen: " << czas_obliczen_NEH.count() << "s" << endl;
    cout << "Calkowity czas obliczen: " << czas_obliczen_NEH.count() << "s" << endl;
    wyjscie.close();
}

long algorytmNEH(int &n, int &m, int **&p, long &cmax, int *&perm) {
    int *priority = new int[n + 1];
    int i, j, l, t; // zmienne pomocnicze
    int **r = new int *[m + 1];
    for (i = 0; i <= m; ++i) r[i] = new int[n + 1]; // dlugosci drogi do wierzch (1,1)
    int **q = new int *[m + 2];
    for (i = 0; i <= m + 1; ++i) q[i] = new int[n + 2]; // dlugosci drogi do wierzch (m,i)
    int *dist = new int[m + 1];
    long c;
    int zad_kryt, temp;

    for (i = 1; i <= n; ++i) { //liczenie priorytetow
        t = 0;
        for (j = 1; j <= m; ++j)
            t += p[j][i];
        priority[i] = t * n - i;  // suma wszystkich zadan
    }
    for (i = 1; i <= n; ++i) // stworz permutacje bazowa, zeby moc sortowac
        perm[i] = i;

    heapSort(priority, perm, n);

    for (i = 0; i <= m; ++i) r[i][0] = 0;  // zadanie 0 nie istnieje
    for (j = 0; j <= n; ++j) r[0][j] = 0;  // maszyna 0 nie istnieje
    for (j = 0; j <= n + 1; ++j) q[m + 1][j] = 0; // maszyna m+1 nie istnieje
    dist[0] = 0; // ze wzgledu na brak definicji

    t = 1;
    for (int zad = 2; zad <= n; ++zad) {
        // odswiez r i q
        for (i = 1; i <= m; ++i)
            for (j = t; j <= zad; ++j)  // do zad, zamiast zad-1, aby policzyć cmax dla kolejnego zadania na koncu
                r[i][j] = max(r[i - 1][j], r[i][j - 1]) + p[i][perm[j]]; // licz odleglosci do (1,1)

        for (i = 0; i <= m; ++i) q[i][zad] = 0; //prawy brzeg 'q' przed dodaniem zadania 'zad'

        for (i = m; i >= 1; --i)
            for (j = zad - 1; j >= t; --j)  // przepisz q z j-1 na j od pozycji t
                q[i][j] = q[i][j - 1];

        for (i = m; i >= 1; --i)
            for (j = t; j >= 1; --j)
                q[i][j] = max(q[i][j + 1], q[i + 1][j]) + p[i][perm[j]];

        // wyznacz najlepsza pozycje w permutacji
        cmax = r[m][zad];                           // cmax najdluzsza sciezka w grafie
        t = zad;                                    // pozycja kolejnego zadania
        for (j = zad - 1; j >= 1; --j) { // wyznacz najlepsza pozycje dla nowego zadania i wylicz cmax
            for (l = 1; l <= m; ++l) // wyznacz D(l), droge z po wstawieniu zadania w j'te miejsce
                dist[l] = max(dist[l - 1], r[l][j - 1]) + p[l][perm[zad]];
            c = dist[1] + q[1][j];
            for (l = 2; l <= m; ++l)
                c = max(c, (long) dist[l] + q[l][j]);
            if (c <= cmax) {
                cmax = c;
                t = j;  // pozycja dla zadania , jesli nie znajdziesz w petli, tzn, ze pozycja na koncu permutacji
            }
        }
        // zapisz permutacje
        i = perm[zad];
        for (j = zad; j > t; --j)
            perm[j] = perm[j - 1];
        perm[t] = i;
        // Modyfikacja IR4
        // po wstawieniu zadania j, znajdz zadanie o najwiekszym wplywie na cmax i znajdz dla niego najlepsza pozycje
        modyfikacjaIR4(m, p, r, q, zad, t, zad_kryt, cmax, perm, dist, temp);
        // Modyfikacja IRR4
        // modyfikacja:wyjecie zadania j'tego (tutaj 'zad') i znalezienie dla niego lepszego miejsca)
        modyfikacjaIRR4(m, p, r, q, zad, t, zad_kryt, cmax, perm, dist, temp);
        t = temp; // zeby zachowac ciaglosc miedzy neh bez modyfikacji i z modyfikacja
    }
    for (i = 0; i <= m + 1; ++i)
        delete[] q[i];
    delete[] q;
    for (i = 0; i <= m; ++i)
        delete[] r[i];
    delete[] r;
    delete[] dist;
    delete[] priority;
    return cmax;
}

void
modyfikacjaIR4(int &m, int **&p, int **&r, int **&q, int &zad, int &t, int &zad_kryt, long &cmax, int *&perm,
               int *&dist, int &temp) {
    int i, j, l;
    int min_cmax, new_cmax;
    long c;
    for (i = 1; i <= m; ++i)
        for (j = t; j <= zad; ++j)
            r[i][j] = max(r[i - 1][j], r[i][j - 1]) + p[i][perm[j]];

    for (i = 0; i <= m; ++i) q[i][zad + 1] = 0;

    for (i = m; i >= 1; --i)
        for (j = zad; j >= t; --j)
            q[i][j] = q[i][j - 1];

    for (i = m; i >= 1; --i)
        for (j = t; j >= 1; --j)
            q[i][j] = max(q[i][j + 1], q[i + 1][j]) + p[i][perm[j]];
    // znajdz najbardziej znaczace zadanie (ktorego usuniecie powoduje najwiekszy spadek cmax)
    for (j = zad, min_cmax = r[m][zad]; j >= 1; --j) {
        if (j != t) {
            new_cmax = 0;
            for (i = m; i >= 1; --i)
                new_cmax = max(new_cmax, r[i][j - 1] + q[i][j + 1]);
            if (new_cmax <= min_cmax) {
                zad_kryt = j;
                min_cmax = new_cmax;
            }
        }
    }
    //odswiez permutacje - zabierz zadanie
    i = perm[zad_kryt];
    for (j = zad_kryt; j <= zad - 1; ++j)
        perm[j] = perm[j + 1];
    perm[zad] = i;
    // odswiez r i q po usunieciu zadania
    for (i = 1; i <= m; ++i)
        for (j = zad_kryt; j <= zad; ++j)
            r[i][j] = max(r[i - 1][j], r[i][j - 1]) + p[i][perm[j]];

    for (i = m; i >= 1; --i) // przepisz q z j-1 na j od pozycji t
        for (j = zad_kryt; j <= zad - 1; ++j)
            q[i][j] = q[i][j + 1];
    for (i = 0; i <= m; ++i) q[i][zad] = 0; // prawy brzeg 'q' przed dodaniem zadania 'zad'

    for (i = m; i >= 1; --i)
        for (j = zad_kryt; j >= 1; --j)
            q[i][j] = max(q[i][j + 1], q[i + 1][j]) + p[i][perm[j]];
    // r i q odswiezone
    // wyznacz najlepsza pozycje dla zadania krytycznego
    cmax = r[m][zad]; // cmax najdluzsza sciezka w grafie
    temp = zad; // miejsce na wyciagniete zadanie
    for (j = zad - 1; j >= 1; --j) {
        for (l = 1; l <= m; ++l)
            dist[l] = max(dist[l - 1], r[l][j - 1]) + p[l][perm[zad]];
        c = dist[1] + q[1][j];
        for (l = 2; l <= m; ++l)
            c = max(c, (long) dist[l] + q[l][j]);
        if (c <= cmax) {
            cmax = c;
            temp = j;  // pozycja dla zadania , jesli nie znajdziesz w petli, tzn, ze pozycja na koncu permutacji
        }
    }
    // odswiez permutacje po modyfikacji pozycji zadania kryt.
    i = perm[zad];
    for (j = zad; j > temp; --j) perm[j] = perm[j - 1];
    perm[temp] = i;

}


void heapSort(int *&priority, int *&perm, int &n) {
    int p_ind;  // rodzic (parent index)
    int ch_ind; //  dziecko (child index)
    int gch_ind; // wieksze dziecko (greater child index)
    int x, y;

    for (int i = 2; i <= n; i++) { // budowanie kopca
        ch_ind = i;
        p_ind = ch_ind / 2;
        x = priority[i];
        y = perm[i];

        while ((p_ind > 0) && (priority[p_ind] > x)) {
            perm[ch_ind] = perm[p_ind];
            priority[ch_ind] = priority[p_ind];
            ch_ind = p_ind;
            p_ind = (ch_ind >> 1);
        }
        priority[ch_ind] = x;
        perm[ch_ind] = y;
    }

    for (int i = n; i > 1; i--)      // rozebranie kopca
    {
        swap(priority[1], priority[i]); // element najwiekszy na koniec
        swap(perm[1], perm[i]);
        p_ind = 1;
        ch_ind = 2;
        while (ch_ind < i) //z
        {
            if ((ch_ind + 1 < i) &&
                (priority[ch_ind + 1] < priority[ch_ind])) //jesli istnieje prawe dziecko i jest mniejsze od lewego
                gch_ind = ch_ind + 1; // mniejsze prawe
            else
                gch_ind = ch_ind;   // mniejsze lewe
            if (priority[gch_ind] >= priority[p_ind])
                break; // jesli mniejsze dziecko jest wieksze od rodzica przerwij 'while'
            swap(priority[p_ind], priority[gch_ind]); // jesli nie, zamien rodzica z wiekszym dzieckiem
            swap(perm[p_ind], perm[gch_ind]);
            p_ind = gch_ind;
            ch_ind = (p_ind << 1);
        }
    } // zmniejsz rozmiar kopca
}

void
modyfikacjaIRR4(int &m, int **&p, int **&r, int **&q, int &zad, int &t, int &zad_kryt, long &cmax, int *&perm,
                int *&dist, int &temp) {
    int i, j, l;
    long c;
// odswiez r i q po wstawieniu zadania krytycznego
    for (i = 1; i <= m; ++i)
        for (j = temp; j <= zad; ++j)
            r[i][j] = max(r[i - 1][j], r[i][j - 1]) + p[i][perm[j]];

    for (i = 0; i <= m; ++i) q[i][zad + 1] = 0; // prawy brzeg 'q' przed dodaniem zadania 'zad'

    for (i = m; i >= 1; --i) // przepisz q z j-1 na j od pozycji t
        for (j = zad; j >= temp; --j)
            q[i][j] = q[i][j - 1];

    for (i = m; i >= 1; --i)
        for (j = temp; j >= 1; --j)
            q[i][j] = max(q[i][j + 1], q[i + 1][j]) + p[i][perm[j]];     // licz odleglosci do(m,zad)

// jesli zmienila sie pozycja zadania t, wyznacz ja
    if (zad_kryt > t && temp <= t) ++t;
    else if (zad_kryt < t && temp >= t) --t;

// odswiez permutacje
    i = perm[t]; // wyciagnij zadanie t z permutacji i umiesc na koncu
    for (j = t; j <= zad - 1; ++j) perm[j] = perm[j + 1];
    perm[zad] = i;

// odswiez r i q po wyciagnieciu zadania
    for (i = 1; i <= m; ++i)
        for (j = t; j <= zad; ++j)
            r[i][j] = max(r[i - 1][j], r[i][j - 1]) + p[i][perm[j]];

    for (i = m; i >= 1; --i) // przepisz q z j-1 na j od pozycji t
        for (j = t; j <= zad - 1; ++j)
            q[i][j] = q[i][j + 1];
    for (i = 0; i <= m; ++i) q[i][zad] = 0; // prawy brzeg 'q' przed dodaniem zadania 'zad'

    for (i = m; i >= 1; --i)
        for (j = t; j >= 1; --j)
            q[i][j] = max(q[i][j + 1], q[i + 1][j]) + p[i][perm[j]];
// r i q odswiezone
// wyznacz najlepsza pozycje dla zadania i oblicz cmax
    cmax = r[m][zad];  // cmax najdluzsza sciezka w grafie
    temp = zad; // nowa pozycja dla zadania
    for (j = zad - 1; j >= 1; --j) {
        for (l = 1; l <= m; ++l)  // wyznacz D(l), droge z po wstawieniu zadania w j'te miejsce
            dist[l] = max(dist[l - 1], r[l][j - 1]) + p[l][perm[zad]];
        c = dist[1] + q[1][j];
        for (l = 2; l <= m; ++l)
            c = max(c, (long) dist[l] + q[l][j]);
        if (c <= cmax) {
            cmax = c;
            temp = j;  // pozycja dla zadania , jesli nie znajdziesz w petli, tzn, ze pozycja na koncu permutacji
        }
    }
// odswiez permutacje
    i = perm[zad];
    for (j = zad; j > temp; --j) perm[j] = perm[j - 1];
    perm[temp] = i;// permutacja po wstawieniu zadanie na pozycji t
}
