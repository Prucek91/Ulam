// g++ lab05_ulam.cpp -pthread -O3 -fopenmp -o binarka
// env OMP_NUM_THREADS=4 ./binarka
#include <iostream>
#include <fstream>
#include <cmath>
#include <omp.h>

using namespace std;



int ** t;
const int size = 1000;
// tablice robocze
const int mx = size;
const int my = size;
const int array_size = size;
// tablica PPM
unsigned char image[size][size][3];
// tablica kolorow
unsigned char colorTheme[][3] = {{255, 100, 255}, {0, 255, 0},{255, 0, 0}, {0, 0, 255}, {154, 206, 233}};

void prepare_ppm()
{
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
        {
            image[i][j][0] = 255;
            image[i][j][1] = 255;
            image[i][j][2] = 255;
        }
}

void save_to_ppm()
{
    FILE *fp;
    char *filename = "ulam.ppm";
    char *comment = "# ";
    fp = fopen(filename, "wb");
    fprintf(fp, "P6\n %s\n %d\n %d\n %d\n", comment, size, size, 255);    
    fwrite(image, 1, 3 * size * size, fp);
    fclose(fp);
    return ;
}


int ulam_get_map(int x, int y, int n) 
{ 
    x -= (n - 1) / 2; 
    y -= n / 2;
    
    int mx = abs(x), my = abs(y); 
    int l = 2 * max(mx, my); 
    int d = y >= x ? l * 3 + x + y : l - x - y; 
    return pow(l - 1, 2) + d; 
} 


int isprime(int n) 
{ 
        int p; 
        for (p = 2; p*p <= n; p++) 
                if (n%p == 0) return 0; 
        return n > 2; 
} 




int main(int ac, char** av)
{
    prepare_ppm();
    t = new int *[mx];
    for (int i = 0; i < mx; i++)
        t[i] = new int[my];

    // ulam     
    int i, j;

    int chunk_size = size / omp_get_max_threads();
    double  start = omp_get_wtime();


    // #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(static)
    // #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(dynamic)
    // #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(guided)
    // #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(static, 25)
    // #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(dynamic, 25)
    // #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(guided, 25)
    #pragma omp parallel for default(shared) shared(t, mx, my) private(i,j) schedule(guided, chunk_size)
    for ( i = 0; i < mx;++i)
    {
        for ( j = 0; j < my; j++)
        {
            t[i][j] = ulam_get_map(i, j, array_size);
            if ( isprime( t[i][j] ) )
            {   
                image[i][j][0] = 0;
                image[i][j][1] = 0;
                image[i][j][2] = 0;
            }
            else
            {
                // get thread ID and color pixels
                int id = omp_get_thread_num();
                image[i][j][0] = colorTheme[id][0];
                image[i][j][1] = colorTheme[id][1];
                image[i][j][2] = colorTheme[id][2];
            }            
        }
    }
    double  stop = omp_get_wtime();
    cout << "Time: " << stop - start << "s"<< endl;
    
    cout << "Numthreads: " << omp_get_max_threads() << endl;

    
   save_to_ppm();
    for (int i = 0; i < my; i++)
        delete [] t[i];
    delete t;
    return 0;
}
