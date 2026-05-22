# Demo Streamlit Genotype Imputation Berbasis HMM ala IMPUTE5

Project ini adalah implementasi edukatif dari genotype imputation berdasarkan
ide utama paper IMPUTE5.

Tujuan project ini **bukan** mereplikasi seluruh software IMPUTE5 asli. Tujuan
utamanya adalah membuat alur algoritma lebih mudah dipahami untuk presentasi:

```text
target haplotype dengan marker hilang
        +
reference haplotype lengkap
        |
        v
seleksi state PBWT-inspired
        |
        v
HMM yang diperkecil
        |
        v
forward-backward posterior probabilities
        |
        v
probabilitas allele hasil imputasi
```

## Cara Menjalankan

Buka PowerShell:

```powershell
cd C:\Tugas\codex_impute5_demo
streamlit run app.py
```

Jika komputer lain belum memiliki package yang dibutuhkan:

```powershell
pip install -r requirements.txt
```

Lalu jalankan lagi:

```powershell
streamlit run app.py
```

## File Dalam Project

- `app.py`: aplikasi Streamlit dan implementasi utama algoritma.
- `simplified_impute5.py`: versi command-line awal dari ide yang sama.
- `requirements.txt`: daftar package untuk menjalankan atau deploy app.

## Apa Yang Ditunjukkan Aplikasi Ini?

Aplikasi ini menunjukkan versi sederhana dari genotype imputation ala IMPUTE5.

### 1. Data Input

Reference panel adalah tabel haplotype lengkap:

```text
h0 = 0 1 0 1 ...
h1 = 0 1 0 1 ...
h2 = 0 1 1 1 ...
```

Target haplotype adalah haplotype yang ingin dilengkapi:

```text
target = 0 1 ? 1 ? 0 ? 1 ? 1
```

Dalam aplikasi:

- `0` dan `1` adalah nilai allele.
- `?` berarti marker pada target hilang dan harus diimputasi.
- kolom `M0`, `M1`, `M2`, ... mewakili marker/SNP.
- baris `h0`, `h1`, `h2`, ... mewakili reference haplotype.

Kamu bisa mengubah reference panel dan target haplotype langsung di aplikasi.

### 2. Mengapa Menggunakan HMM?

Dalam genotype imputation, target haplotype dapat dianggap sebagai gabungan atau
mosaic dari reference haplotype yang sudah diketahui.

Dalam istilah HMM:

- hidden state = reference haplotype yang sedang disalin oleh target
- emission probability = seberapa cocok allele target dengan allele reference
- transition probability = peluang berpindah dari satu reference haplotype ke haplotype lain
- posterior probability = probabilitas akhir setiap hidden state setelah forward-backward

### 3. Mengapa HMM Standar Menjadi Berat?

HMM standar memakai semua reference haplotype sebagai hidden state.

Contoh:

```text
10.000 reference haplotypes = 10.000 HMM states
1.000.000 reference haplotypes = 1.000.000 HMM states
```

Semakin banyak state, semakin berat komputasi forward-backward.

### 4. Ide Utama IMPUTE5

IMPUTE5 tetap memakai HMM, tetapi membuatnya lebih scalable.

Daripada memakai semua reference haplotype sebagai state, IMPUTE5 memakai PBWT
untuk mencari reference haplotype yang mirip secara lokal dengan target. HMM
kemudian hanya dijalankan pada subset yang terpilih.

Dalam demo ini, PBWT asli disederhanakan menjadi seleksi neighbour
PBWT-inspired:

1. Lihat marker target yang sudah teramati.
2. Urutkan reference haplotype berdasarkan pola allele terdekat.
3. Masukkan pola target ke urutan tersebut.
4. Pilih reference haplotype terdekat sebagai kandidat state.

Intinya:

```text
HMM standar:
semua reference haplotype menjadi state

HMM ala IMPUTE5:
hanya haplotype yang mirip secara lokal menjadi state
```

### 5. Gambaran PBWT Asli

Demo Streamlit memakai cara sederhana: menghitung kecocokan allele pada marker
target yang sudah diketahui.

PBWT asli di IMPUTE5 bekerja lebih canggih:

1. Pada setiap marker, reference haplotype diurutkan berdasarkan pola allele
   sebelumnya, dibaca dari marker saat ini ke arah kiri. Ini disebut reverse prefix.
2. Haplotype yang berdekatan dalam urutan PBWT biasanya memiliki segmen allele
   yang panjang dan mirip.
3. Target haplotype disisipkan secara konseptual ke urutan tersebut.
4. Reference haplotype yang dekat dengan posisi target dipilih sebagai kandidat
   copying state.

Contoh sederhana:

```text
target pada marker teramati:
M0=0, M1=1, M3=1, M5=0

Jika sedang berada di M5, PBWT melihat pola dari kanan ke kiri:
M5, M3, M1, M0 = 0, 1, 1, 0
```

Reference yang memiliki pola kanan-ke-kiri mirip akan berada dekat dengan target
dalam urutan PBWT. Itulah alasan reference tersebut dipilih.

Perbedaan ringkas:

```text
Demo:
hitung kecocokan sederhana agar mudah dipahami

IMPUTE5 asli:
pakai PBWT/FM-index untuk mencari neighbour secara cepat pada data sangat besar
```

### 6. HMM Asli Dalam IMPUTE5

IMPUTE5 memakai model probabilistik keluarga Li-and-Stephens HMM.

Hidden state:

```text
state = reference haplotype yang sedang disalin oleh target
```

Emission probability:

```text
seberapa masuk akal allele target jika target sedang menyalin dari reference tertentu
```

Jika allele target cocok dengan allele reference, probabilitasnya tinggi. Jika
tidak cocok, probabilitasnya rendah, tetapi tidak nol.

Transition probability:

```text
peluang target tetap menyalin reference yang sama
atau berpindah ke reference haplotype lain
```

Pada data asli, peluang berpindah dipengaruhi genetic map dan jarak antar marker:

- marker dekat -> peluang tetap pada state yang sama lebih tinggi
- marker jauh -> peluang berpindah state lebih besar

Forward-backward:

```text
menggabungkan informasi dari kiri dan kanan marker
untuk menghitung posterior probability setiap state
```

Setelah posterior didapat:

```text
P(allele = 1) = jumlah posterior state yang reference allele-nya 1
```

Jadi angka probabilitas imputasi berasal dari:

```text
kecocokan allele target-reference
+ peluang switch antar haplotype
+ informasi marker kiri dan kanan
+ allele reference pada marker yang hilang
```

### 7. Forward-Backward

Setelah state dipilih, aplikasi menjalankan forward-backward:

- forward pass memakai informasi dari sisi kiri sequence
- backward pass memakai informasi dari sisi kanan sequence
- keduanya digabung menjadi posterior copying probabilities

Untuk setiap marker yang hilang, aplikasi mengecek reference haplotype terpilih
mana yang memiliki allele `1` pada marker tersebut. Probabilitas state-state itu
dijumlahkan untuk menghasilkan:

```text
P(allele = 1)
```

Jadi output bukan hanya:

```text
marker = 0 atau 1
```

Tetapi:

```text
P(allele = 1) = 0.683
prediksi allele = 1
```

Ini penting karena genotype imputation bersifat probabilistik.

## Kontrol Interaktif

Aplikasi Streamlit memungkinkan kamu mengubah:

- nilai reference haplotype
- nilai target haplotype
- marker target yang hilang menggunakan `?`
- jumlah neighbour yang dipilih per marker teramati
- emission error rate
- recombination/switch rate

Nilai default sudah dibuat sebagai contoh yang langsung bisa berjalan.

## Bagaimana Kalau Memakai IMPUTE5 Asli?

IMPUTE5 asli dapat dipakai, tetapi sifatnya berbeda dari demo ini.

Demo ini:

```text
Python/Streamlit
data kecil 0/1/?
tujuan: menjelaskan algoritma
```

IMPUTE5 asli:

```text
software command-line
data genomik nyata
tujuan: pipeline bioinformatika skala besar
```

Pipeline IMPUTE5 asli secara umum:

```text
phased target genotype file
+ reference panel
+ genetic map
+ region kromosom
        |
        v
impute5 command-line tool
        |
        v
output genotype hasil imputasi
```

Contoh bentuk command secara umum:

```bash
impute5 \
  --h reference_panel.bcf \
  --g target_panel.bcf \
  --m genetic_map.txt \
  --r chr20:1000000-2000000 \
  --o imputed_output.bgen
```

Catatan: opsi persisnya bisa berbeda tergantung versi IMPUTE5 dan format input
yang digunakan.

Jika memakai IMPUTE5 asli, biasanya dibutuhkan:

- software/binary IMPUTE5,
- Linux, WSL, atau server/HPC,
- target genotype yang sudah di-phase,
- reference panel,
- genetic map,
- file genomik seperti VCF, BCF, BGEN, atau imp5,
- pemahaman format kromosom dan region.

Untuk presentasi implementasi, demo Streamlit lebih cocok karena algoritmanya
terlihat dan bisa dijelaskan. IMPUTE5 asli lebih cocok untuk menjalankan analisis
genomik nyata.

## Penyederhanaan Dibanding IMPUTE5 Asli

Demo ini sengaja dibuat kecil dan mudah dibaca.

IMPUTE5 asli memiliki:

- reference panel sangat besar,
- data phased diploid,
- PBWT/FM-index yang dioptimalkan,
- genomic window,
- genetic map,
- format BGEN dan `imp5`,
- optimasi multiprocessing dan multithreading.

Demo ini mempertahankan ide utama:

```text
pilih state HMM yang relevan dulu sebelum forward-backward dijalankan
```

## Ringkasan

> HMM standar memakai semua reference haplotype sebagai hidden state. IMPUTE5
> mempercepat proses dengan memakai PBWT untuk memilih haplotype yang mirip
> secara lokal terlebih dahulu. Setelah itu, HMM dijalankan pada state yang lebih
> sedikit, sehingga proses imputasi menjadi lebih efisien tanpa menghilangkan
> struktur probabilistik dari model Li-and-Stephens.
