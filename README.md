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
        |
        v
haplotype lengkap (sebelum vs sesudah imputasi)
```

## Cara Menjalankan

Buka terminal (PowerShell atau Codespace):

```bash
pip install -r requirements.txt
streamlit run app.py
```

Jika `streamlit` tidak ditemukan setelah install, coba:

```bash
python -m streamlit run app.py
```

Untuk GitHub Codespace, setelah Streamlit berjalan:
- Buka tab **Ports** di bagian bawah
- Cari port **8501** → klik kanan → **Port Visibility → Public**
- Buka URL Codespace (`-8501.app.github.dev`) di tab browser baru secara manual

## File Dalam Project

- `app.py`: aplikasi Streamlit dan implementasi utama algoritma.
- `requirements.txt`: daftar package untuk menjalankan atau deploy app.

## Apa Yang Ditunjukkan Aplikasi Ini?

Aplikasi ini menunjukkan versi sederhana dari genotype imputation ala IMPUTE5,
dibagi menjadi empat bagian utama.

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
- Kolom `M0`, `M1`, `M2`, ... mewakili marker/SNP.
- Baris `h0`, `h1`, `h2`, ... mewakili reference haplotype.

Kamu bisa mengubah reference panel dan target haplotype langsung di aplikasi.
Semua penjelasan dan contoh angka di bawahnya akan ikut berubah secara otomatis.

### 2. Seleksi State PBWT-Inspired

Sebelum HMM dijalankan, aplikasi mencari reference haplotype yang paling mirip
dengan target pada marker yang sudah diketahui.

Langkah-langkahnya:

1. Ambil pola target hanya pada marker teramati (bukan `?`).
2. Bandingkan pola tiap reference terhadap pola target — hitung berapa posisi yang cocok.
3. Urutkan reference dari yang paling mirip.
4. Ambil sejumlah reference teratas (diatur oleh slider) sebagai kandidat state HMM.

Expander **"Apa maksud marker teramati dan pola allele terdekat?"** menjelaskan
proses ini secara langkah demi langkah dengan contoh angka dari data yang sedang aktif.

### 3. HMM Yang Diperkecil

Setelah state terpilih, HMM dijalankan hanya pada subset reference tersebut.

- **Emission probability** — seberapa cocok allele target dengan allele reference di marker ini.
  Jika cocok, probabilitas tinggi (`1 - error_rate`). Jika tidak cocok, probabilitas rendah (`error_rate`).
- **Transition probability** — peluang model tetap menyalin dari reference yang sama atau
  berpindah ke reference lain (merepresentasikan rekombinasi).
- **Forward-backward** — menggabungkan informasi dari sisi kiri dan kanan setiap marker
  untuk menghasilkan posterior probability tiap state.

Expander **"Dari mana angka emission, transition, dan posterior muncul?"** menjelaskan
ketiga konsep ini dengan angka konkret yang dihitung dari nilai slider dan data saat ini.

### 4. Hasil Imputasi

Untuk setiap marker yang hilang (`?`), aplikasi menghitung:

```text
P(allele = 1) = jumlah posterior state yang reference allele-nya 1
```

Output akhir mencakup:

- Tabel probabilitas per marker yang hilang.
- Bar chart `P(allele = 1)`.
- Contoh cara membaca angka untuk marker pertama yang hilang.
- **Tabel perbandingan haplotype sebelum dan sesudah imputasi** — menunjukkan
  haplotype target lengkap dengan semua `?` sudah diganti prediksi akhir.
- Tabel ringkasan semua marker dengan status "sudah diketahui" atau "diimputasi".

## Kontrol Interaktif

Aplikasi Streamlit memungkinkan kamu mengubah:

- Nilai reference haplotype dan target haplotype (langsung di tabel)
- Marker target yang hilang menggunakan `?`
- Jumlah neighbour PBWT-inspired per marker teramati
- Emission error rate
- Recombination / switch rate

Semua penjelasan teks, contoh angka, dan tabel akan otomatis menyesuaikan
dengan perubahan input.

## Mengapa Menggunakan HMM?

Dalam genotype imputation, target haplotype dapat dianggap sebagai gabungan atau
mosaic dari reference haplotype yang sudah diketahui.

Dalam istilah HMM:

- **Hidden state** = reference haplotype yang sedang disalin oleh target
- **Emission probability** = seberapa cocok allele target dengan allele reference
- **Transition probability** = peluang berpindah dari satu reference haplotype ke haplotype lain
- **Posterior probability** = probabilitas akhir setiap hidden state setelah forward-backward

### Mengapa HMM Standar Menjadi Berat?

```text
10.000 reference haplotypes  = 10.000 HMM states
1.000.000 reference haplotypes = 1.000.000 HMM states
```

### Ide Utama IMPUTE5

Daripada memakai semua reference haplotype sebagai state, IMPUTE5 memakai PBWT
untuk mencari reference haplotype yang mirip secara lokal dengan target terlebih
dahulu. HMM kemudian hanya dijalankan pada subset yang terpilih.

```text
HMM standar:    semua reference haplotype menjadi state
HMM ala IMPUTE5: hanya haplotype yang mirip secara lokal menjadi state
```

## Penyederhanaan Dibanding IMPUTE5 Asli

Demo ini sengaja dibuat kecil dan mudah dibaca.

| Aspek | Demo Streamlit | IMPUTE5 asli |
|---|---|---|
| Bahasa/platform | Python/Streamlit | Software command-line (Linux/WSL/HPC) |
| Data | Array kecil `0/1/?` | VCF, BCF, BGEN, imp5 |
| Seleksi state | Hitung kecocokan sederhana | PBWT/FM-index yang dioptimalkan |
| Skala | Beberapa haplotype | Ribuan hingga jutaan haplotype |
| Tujuan | Menjelaskan algoritma | Pipeline bioinformatika skala besar |

## Bagaimana Kalau Memakai IMPUTE5 Asli?

Pipeline IMPUTE5 asli secara umum:

```bash
impute5 \
  --h reference_panel.bcf \
  --g target_panel.bcf \
  --m genetic_map.txt \
  --r chr20:1000000-2000000 \
  --o imputed_output.bgen
```

Yang dibutuhkan: binary IMPUTE5, Linux/WSL/HPC, data phased diploid,
reference panel, genetic map, dan format file genomik (VCF/BCF/BGEN/imp5).

Untuk presentasi implementasi, demo Streamlit lebih cocok. IMPUTE5 asli
lebih cocok untuk analisis genomik nyata.

## Ringkasan

> HMM standar memakai semua reference haplotype sebagai hidden state. IMPUTE5
> mempercepat proses dengan memakai PBWT untuk memilih haplotype yang mirip
> secara lokal terlebih dahulu. Setelah itu, HMM dijalankan pada state yang lebih
> sedikit, sehingga proses imputasi menjadi lebih efisien tanpa menghilangkan
> struktur probabilistik dari model Li-and-Stephens.