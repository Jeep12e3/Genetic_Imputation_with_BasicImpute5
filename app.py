from bisect import bisect_left
from typing import Iterable, Optional

import pandas as pd
import streamlit as st


Allele = int
MaybeAllele = Optional[int]


DEFAULT_REFERENCE_PANEL: list[list[Allele]] = [
    [0, 1, 0, 1, 1, 0, 0, 1, 0, 1],
    [0, 1, 0, 1, 0, 0, 0, 1, 0, 1],
    [0, 1, 1, 1, 1, 0, 1, 1, 0, 1],
    [1, 0, 1, 0, 0, 1, 1, 0, 1, 0],
    [1, 0, 1, 0, 1, 1, 1, 0, 1, 0],
    [0, 1, 0, 1, 1, 0, 0, 0, 0, 1],
    [1, 1, 0, 1, 1, 0, 0, 1, 1, 1],
    [0, 0, 1, 0, 0, 1, 1, 0, 0, 0],
]

DEFAULT_TARGET: list[MaybeAllele] = [0, 1, None, 1, None, 0, None, 1, None, 1]


def marker_columns(marker_count: int) -> list[str]:
    return [f"M{index}" for index in range(marker_count)]


def default_reference_dataframe() -> pd.DataFrame:
    return pd.DataFrame(
        DEFAULT_REFERENCE_PANEL,
        columns=marker_columns(len(DEFAULT_REFERENCE_PANEL[0])),
        index=[f"h{index}" for index in range(len(DEFAULT_REFERENCE_PANEL))],
    )


def default_target_dataframe() -> pd.DataFrame:
    values = ["?" if allele is None else str(allele) for allele in DEFAULT_TARGET]
    return pd.DataFrame(
        [values],
        columns=marker_columns(len(DEFAULT_TARGET)),
        index=["target"],
    )


def parse_reference_dataframe(dataframe: pd.DataFrame) -> list[list[Allele]]:
    return dataframe.astype(int).values.tolist()


def parse_target_dataframe(dataframe: pd.DataFrame) -> list[MaybeAllele]:
    target: list[MaybeAllele] = []
    for value in dataframe.iloc[0].tolist():
        if str(value) == "?":
            target.append(None)
        else:
            target.append(int(value))
    return target


def reverse_prefix_key(
    haplotype: list[MaybeAllele],
    observed_positions: list[int],
    marker: int,
) -> tuple[int, ...]:
    prefix_positions = [position for position in observed_positions if position <= marker]
    return tuple(int(haplotype[position]) for position in reversed(prefix_positions))


def choose_neighbours(
    sorted_items: list[tuple[tuple[int, ...], int]],
    insert_at: int,
    count: int,
) -> list[int]:
    chosen: list[int] = []
    left = insert_at - 1
    right = insert_at

    while len(chosen) < count and (left >= 0 or right < len(sorted_items)):
        if left >= 0:
            chosen.append(sorted_items[left][1])
            left -= 1
            if len(chosen) == count:
                break

        if right < len(sorted_items):
            chosen.append(sorted_items[right][1])
            right += 1

    return chosen


def pbwt_inspired_state_selection(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    neighbours_per_marker: int,
) -> list[int]:
    observed_positions = [index for index, allele in enumerate(target) if allele is not None]
    if not observed_positions:
        return list(range(len(reference_panel)))

    selected: list[int] = []

    for marker in observed_positions:
        target_key = reverse_prefix_key(target, observed_positions, marker)
        sorted_items = sorted(
            (
                reverse_prefix_key(reference_haplotype, observed_positions, marker),
                reference_index,
            )
            for reference_index, reference_haplotype in enumerate(reference_panel)
        )
        keys_only = [item[0] for item in sorted_items]
        insert_at = bisect_left(keys_only, target_key)
        selected.extend(choose_neighbours(sorted_items, insert_at, neighbours_per_marker))

    unique_selected: list[int] = []
    seen: set[int] = set()
    for state in selected:
        if state not in seen:
            unique_selected.append(state)
            seen.add(state)

    return unique_selected


def emission_probability(
    observed_allele: MaybeAllele,
    reference_allele: Allele,
    error_rate: float,
) -> float:
    if observed_allele is None:
        return 1.0
    return 1.0 - error_rate if observed_allele == reference_allele else error_rate


def transition_probability(
    previous_state: int,
    current_state: int,
    number_of_states: int,
    recombination_rate: float,
) -> float:
    if number_of_states == 1:
        return 1.0
    if previous_state == current_state:
        return 1.0 - recombination_rate
    return recombination_rate / (number_of_states - 1)


def normalize(values: Iterable[float]) -> list[float]:
    values = list(values)
    total = sum(values)
    if total == 0.0:
        return [1.0 / len(values)] * len(values)
    return [value / total for value in values]


def forward_algorithm(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    error_rate: float,
    recombination_rate: float,
) -> list[list[float]]:
    marker_count = len(target)
    state_count = len(selected_states)
    alpha = [[0.0] * state_count for _ in range(marker_count)]

    alpha[0] = normalize(
        (1.0 / state_count)
        * emission_probability(target[0], reference_panel[ref_index][0], error_rate)
        for ref_index in selected_states
    )

    for marker in range(1, marker_count):
        row = []
        for current_state_position, ref_index in enumerate(selected_states):
            incoming_probability = sum(
                alpha[marker - 1][previous_state_position]
                * transition_probability(
                    previous_state_position,
                    current_state_position,
                    state_count,
                    recombination_rate,
                )
                for previous_state_position in range(state_count)
            )
            row.append(
                incoming_probability
                * emission_probability(target[marker], reference_panel[ref_index][marker], error_rate)
            )
        alpha[marker] = normalize(row)

    return alpha


def backward_algorithm(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    error_rate: float,
    recombination_rate: float,
) -> list[list[float]]:
    marker_count = len(target)
    state_count = len(selected_states)
    beta = [[0.0] * state_count for _ in range(marker_count)]
    beta[-1] = [1.0] * state_count

    for marker in range(marker_count - 2, -1, -1):
        row = []
        for current_state_position in range(state_count):
            probability = sum(
                transition_probability(
                    current_state_position,
                    next_state_position,
                    state_count,
                    recombination_rate,
                )
                * emission_probability(
                    target[marker + 1],
                    reference_panel[next_ref_index][marker + 1],
                    error_rate,
                )
                * beta[marker + 1][next_state_position]
                for next_state_position, next_ref_index in enumerate(selected_states)
            )
            row.append(probability)
        beta[marker] = normalize(row)

    return beta


def posterior_probabilities(alpha: list[list[float]], beta: list[list[float]]) -> list[list[float]]:
    return [
        normalize(
            alpha[marker][state] * beta[marker][state]
            for state in range(len(alpha[marker]))
        )
        for marker in range(len(alpha))
    ]


def impute_missing_alleles(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    posterior: list[list[float]],
) -> list[Optional[float]]:
    imputed: list[Optional[float]] = []

    for marker, observed_allele in enumerate(target):
        if observed_allele is not None:
            imputed.append(None)
            continue

        probability_allele_1 = sum(
            posterior[marker][state_position]
            for state_position, ref_index in enumerate(selected_states)
            if reference_panel[ref_index][marker] == 1
        )
        imputed.append(probability_allele_1)

    return imputed


def run_imputation(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    neighbours_per_marker: int,
    error_rate: float,
    recombination_rate: float,
) -> tuple[list[int], list[list[float]], list[Optional[float]]]:
    selected_states = pbwt_inspired_state_selection(
        reference_panel,
        target,
        neighbours_per_marker,
    )
    alpha = forward_algorithm(
        reference_panel,
        target,
        selected_states,
        error_rate,
        recombination_rate,
    )
    beta = backward_algorithm(
        reference_panel,
        target,
        selected_states,
        error_rate,
        recombination_rate,
    )
    posterior = posterior_probabilities(alpha, beta)
    imputed = impute_missing_alleles(reference_panel, target, selected_states, posterior)
    return selected_states, posterior, imputed


def reference_display_dataframe(
    reference_panel: list[list[Allele]],
    selected_states: list[int],
) -> pd.DataFrame:
    dataframe = pd.DataFrame(
        reference_panel,
        columns=marker_columns(len(reference_panel[0])),
        index=[f"h{index}" for index in range(len(reference_panel))],
    )
    dataframe.insert(
        0,
        "Dipakai HMM?",
        ["ya" if index in selected_states else "tidak" for index in range(len(reference_panel))],
    )
    return dataframe


def target_display_dataframe(target: list[MaybeAllele]) -> pd.DataFrame:
    return pd.DataFrame(
        [["?" if allele is None else allele for allele in target]],
        columns=marker_columns(len(target)),
        index=["target"],
    )


def result_dataframe(target: list[MaybeAllele], imputed: list[Optional[float]]) -> pd.DataFrame:
    rows = []
    for marker, probability in enumerate(imputed):
        if probability is None:
            continue
        rows.append(
            {
                "Marker": f"M{marker}",
                "Keterangan": "hilang pada target",
                "P(allele = 1)": round(probability, 3),
                "P(allele = 0)": round(1.0 - probability, 3),
                "Prediksi allele": 1 if probability >= 0.5 else 0,
            }
        )
    return pd.DataFrame(rows)


def explain_first_missing_marker(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    posterior: list[list[float]],
    imputed: list[Optional[float]],
) -> Optional[tuple[int, pd.DataFrame, float]]:
    for marker, probability in enumerate(imputed):
        if probability is None:
            continue

        rows = []
        for state_position, reference_index in enumerate(selected_states):
            allele = reference_panel[reference_index][marker]
            posterior_value = posterior[marker][state_position]
            rows.append(
                {
                    "Reference": f"h{reference_index}",
                    f"Allele di M{marker}": allele,
                    "Bobot/posterior": round(posterior_value, 3),
                    "Masuk penjumlahan P(allele=1)?": "ya" if allele == 1 else "tidak",
                }
            )

        return marker, pd.DataFrame(rows), float(probability)

    return None


def observed_marker_dataframe(target: list[MaybeAllele]) -> pd.DataFrame:
    rows = []
    for marker, allele in enumerate(target):
        rows.append(
            {
                "Marker": f"M{marker}",
                "Nilai target": "?" if allele is None else allele,
                "Dipakai untuk mencari kemiripan?": "ya" if allele is not None else "tidak",
                "Alasan": "marker sudah diketahui" if allele is not None else "marker masih hilang",
            }
        )
    return pd.DataFrame(rows)


def similarity_dataframe(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
) -> pd.DataFrame:
    observed_positions = [index for index, allele in enumerate(target) if allele is not None]
    rows = []

    for reference_index, haplotype in enumerate(reference_panel):
        matches = sum(1 for marker in observed_positions if haplotype[marker] == target[marker])
        pattern = " ".join(str(haplotype[marker]) for marker in observed_positions)
        target_pattern = " ".join(str(target[marker]) for marker in observed_positions)
        rows.append(
            {
                "Reference": f"h{reference_index}",
                "Pola pada marker teramati": pattern,
                "Pola target": target_pattern,
                "Jumlah cocok": matches,
            }
        )

    return pd.DataFrame(rows).sort_values("Jumlah cocok", ascending=False)


def probability_example_dataframe(
    reference_panel: list[list[Allele]],
    target: list[MaybeAllele],
    selected_states: list[int],
    marker: int,
    error_rate: float,
    recombination_rate: float,
) -> pd.DataFrame:
    rows = []
    state_count = len(selected_states)

    for current_position, reference_index in enumerate(selected_states):
        reference_allele = reference_panel[reference_index][marker]
        observed_allele = target[marker]
        emission = emission_probability(observed_allele, reference_allele, error_rate)
        stay_probability = transition_probability(
            current_position,
            current_position,
            state_count,
            recombination_rate,
        )
        switch_probability = (
            transition_probability(
                current_position,
                0 if current_position != 0 else min(1, state_count - 1),
                state_count,
                recombination_rate,
            )
            if state_count > 1
            else 0
        )
        rows.append(
            {
                "State": f"h{reference_index}",
                f"Allele reference di M{marker}": reference_allele,
                "Allele target": "?" if observed_allele is None else observed_allele,
                "Emission": round(emission, 3),
                "P(tetap di state ini)": round(stay_probability, 3),
                "P(pindah ke state lain)": round(switch_probability, 3),
            }
        )

    return pd.DataFrame(rows)


def posterior_dataframe(
    selected_states: list[int],
    posterior: list[list[float]],
) -> pd.DataFrame:
    dataframe = pd.DataFrame(
        posterior,
        columns=[f"h{state}" for state in selected_states],
        index=marker_columns(len(posterior)),
    )
    return dataframe.round(3)


def reset_data() -> None:
    st.session_state.reference_editor = default_reference_dataframe()
    st.session_state.target_editor = default_target_dataframe()


st.set_page_config(
    page_title="Demo HMM IMPUTE5",
    page_icon="DNA",
    layout="wide",
)

if "reference_editor" not in st.session_state:
    st.session_state.reference_editor = default_reference_dataframe()
if "target_editor" not in st.session_state:
    st.session_state.target_editor = default_target_dataframe()

st.title("Demo Genotype Imputation Berbasis HMM ala IMPUTE5")
st.caption("Aplikasi Python/Streamlit untuk menjelaskan ide implementasi utama dari paper IMPUTE5.")

st.markdown(
    """
    **Masalah.** Dalam banyak studi genetika, sampel target tidak memiliki semua
    marker/SNP yang mungkin diamati. Genotype imputation digunakan untuk
    memperkirakan marker yang hilang dengan membandingkan target haplotype
    terhadap reference panel yang lebih lengkap.

    **Ide utama.** HMM biasa dapat memakai semua reference haplotype sebagai
    hidden state. Masalahnya, jumlah state menjadi sangat besar jika reference
    panel berisi ribuan sampai jutaan haplotype. IMPUTE5 mengurangi beban ini
    dengan memilih haplotype yang paling mirip secara lokal terlebih dahulu,
    lalu menjalankan HMM hanya pada subset yang terpilih.
    """
)

st.info(
    "Cara membaca demo ini secara sederhana: angka `0` dan `1` adalah dua kemungkinan nilai allele. "
    "Tanda `?` berarti nilai pada marker itu belum diketahui. Aplikasi ini mencoba menebak `?` "
    "dengan melihat baris reference yang polanya mirip dengan target."
)

with st.sidebar:
    st.header("Kontrol")
    st.write("Nilai-nilai ini bisa diubah untuk melihat bagaimana algoritma bereaksi.")
    neighbours_per_marker = st.slider(
        "Jumlah neighbour PBWT-inspired per marker teramati",
        1,
        6,
        4,
        help=(
            "Semakin besar nilainya, semakin banyak reference haplotype yang dipilih "
            "sebelum HMM dijalankan. State lebih banyak bisa memberi cakupan lebih luas, "
            "tetapi komputasi juga lebih berat."
        ),
    )
    error_rate = st.slider(
        "Emission error rate",
        0.001,
        0.100,
        0.010,
        step=0.001,
        format="%.3f",
        help=(
            "Probabilitas bahwa allele target yang diamati berbeda dari allele reference "
            "yang sedang disalin. Dalam demo ini, ini mewakili noise atau ketidakcocokan model."
        ),
    )
    recombination_rate = st.slider(
        "Recombination / switch rate",
        0.01,
        0.30,
        0.08,
        step=0.01,
        help=(
            "Probabilitas berpindah dari satu reference haplotype ke reference haplotype lain. "
            "Dalam konteks genetika, perpindahan ini mewakili rekombinasi sepanjang kromosom."
        ),
    )
    st.button("Reset data default", on_click=reset_data)

with st.expander("Kosakata penting dalam demo ini", expanded=True):
    st.markdown(
        """
        - **Allele**: nilai pada satu marker genetik. Demo sederhana ini hanya memakai `0` dan `1`.
        - **Marker / SNP**: satu posisi pada genom. Kolom `M0`, `M1`, ... mewakili marker.
        - **Haplotype**: urutan allele pada beberapa marker. Baris `h0`, `h1`, ... adalah reference haplotype.
        - **Reference panel**: kumpulan haplotype lengkap yang dipakai sebagai pembanding.
        - **Target haplotype**: haplotype yang ingin kita lengkapi. Simbol `?` berarti marker tersebut hilang.
        - **Hidden state**: dalam HMM, state berarti "target sedang menyalin dari reference haplotype ini".
        - **Posterior probability**: probabilitas akhir bahwa suatu state menjelaskan marker tertentu setelah forward-backward.
        """
    )

with st.expander("Penjelasan sederhana: angka-angka ini datang dari mana?", expanded=True):
    st.markdown(
        """
        Anggap setiap haplotype seperti **baris jawaban** berisi angka `0` dan `1`.
        Reference panel adalah kumpulan baris yang sudah lengkap. Target adalah baris
        yang sebagian angkanya hilang.

        Contoh:

        ```text
        h0     = 0 1 0 1 1 0 0 1 0 1
        h1     = 0 1 0 1 0 0 0 1 0 1
        target = 0 1 ? 1 ? 0 ? 1 ? 1
        ```

        Artinya, target sudah diketahui pada beberapa posisi, tetapi masih memiliki
        nilai `?` pada marker tertentu. Tugas algoritma adalah memperkirakan apakah
        setiap `?` lebih mungkin bernilai `0` atau `1`.

        **Penting:** data `0/1/?` di demo ini adalah data contoh buatan, bukan data
        asli dari paper. Data dibuat kecil supaya cara kerja algoritma mudah dilihat.
        Pada data genomik asli, tabelnya bisa berisi ribuan sampai jutaan haplotype.
        """
    )

st.subheader("1. Ubah Data Input")
st.write(
    "Kamu bisa mengubah reference panel dan target haplotype di bawah ini. "
    "Reference hanya boleh berisi `0` atau `1`. Untuk target, gunakan `?` pada marker yang ingin diimputasi."
)

reference_editor = st.data_editor(
    st.session_state.reference_editor,
    column_config={
        column: st.column_config.SelectboxColumn(column, options=[0, 1], required=True)
        for column in st.session_state.reference_editor.columns
    },
    width="stretch",
    key="reference_data_editor",
)

target_editor = st.data_editor(
    st.session_state.target_editor,
    column_config={
        column: st.column_config.SelectboxColumn(column, options=["0", "1", "?"], required=True)
        for column in st.session_state.target_editor.columns
    },
    width="stretch",
    key="target_data_editor",
)

reference_panel = parse_reference_dataframe(reference_editor)
target = parse_target_dataframe(target_editor)

selected_states, posterior, imputed = run_imputation(
    reference_panel,
    target,
    neighbours_per_marker,
    error_rate,
    recombination_rate,
)

metric_left, metric_mid, metric_right = st.columns(3)
metric_left.metric("Jumlah reference haplotype", len(reference_panel))
metric_mid.metric("State pada HMM standar", len(reference_panel))
metric_right.metric("State terpilih pada demo", len(selected_states))

st.subheader("2. Seleksi State PBWT-Inspired")
st.write(
    "Versi sederhananya: sebelum menebak nilai yang hilang, aplikasi mencari dulu baris reference "
    "yang bentuknya paling mirip dengan target pada bagian yang sudah diketahui. Baris yang "
    "mirip itulah yang dipakai oleh HMM."
)
st.markdown(
    """
    IMPUTE5 asli memakai PBWT untuk menemukan reference haplotype yang memiliki
    kecocokan lokal panjang dengan target. Demo ini memakai pendekatan yang lebih
    mudah dibaca:

    1. Lihat hanya marker target yang teramati, bukan `?`.
    2. Pada setiap marker teramati, urutkan reference haplotype berdasarkan pola allele terdekat.
    3. Tempatkan pola target ke dalam urutan tersebut.
    4. Ambil reference haplotype yang posisinya dekat sebagai kandidat copying state.

    Setelah itu, HMM hanya dijalankan pada state yang terpilih.
    """
)
st.write("State yang terpilih: " + ", ".join(f"h{state}" for state in selected_states))
st.dataframe(reference_display_dataframe(reference_panel, selected_states), width="stretch")

with st.expander("Apa maksud marker teramati dan pola allele terdekat?", expanded=True):
    _observed_pos = [i for i, a in enumerate(target) if a is not None]
    _missing_pos = [i for i, a in enumerate(target) if a is None]
    _target_display = "  ".join("?" if a is None else str(a) for a in target)
    _marker_label = "  ".join(f"M{i}" for i in range(len(target)))
    _observed_label = ", ".join(f"`M{i}={target[i]}`" for i in _observed_pos)
    _missing_label = ", ".join(f"`M{i}`" for i in _missing_pos)

    _obs_header = "  ".join(f"M{i}" for i in _observed_pos)
    _obs_values = "  ".join(str(target[i]) for i in _observed_pos)

    # Build comparison lines for top 2 best and 1 worst reference
    _sim_df = similarity_dataframe(reference_panel, target)
    _total_obs = len(_observed_pos)
    _example_lines = []
    for _, row in _sim_df.head(2).iterrows():
        ref_name = row["Reference"]
        ref_idx = int(ref_name[1:])
        pattern = row["Pola pada marker teramati"]
        matches = row["Jumlah cocok"]
        _example_lines.append(f"{ref_name} pada marker teramati: {pattern}  → {matches}/{_total_obs} cocok")
    _worst_row = _sim_df.iloc[-1]
    _worst_name = _worst_row["Reference"]
    _worst_pattern = _worst_row["Pola pada marker teramati"]
    _worst_matches = _worst_row["Jumlah cocok"]
    _example_lines.append(f"{_worst_name} pada marker teramati: {_worst_pattern}  → {_worst_matches}/{_total_obs} cocok  (paling berbeda)")
    _example_block = "\n        ".join(_example_lines)

    st.markdown(
        f"""
        ### Apa itu marker teramati?

        **Marker teramati** adalah marker pada target yang nilainya sudah diketahui, yaitu bukan `?`.

        Target saat ini:

        ```text
        target = {_target_display}
                 {_marker_label}
        ```

        Marker yang **teramati**: {_observed_label}

        Marker yang **hilang** (tidak dipakai untuk mencari kemiripan): {_missing_label}

        Kenapa yang hilang tidak dipakai? Karena kita belum tahu nilainya — justru itulah yang
        ingin kita tebak. Kita hanya bisa membandingkan apa yang sudah diketahui.
        """
    )
    st.dataframe(observed_marker_dataframe(target), width="stretch", hide_index=True)
    st.markdown(
        f"""
        ---
        ### Bagaimana sistem menentukan "pola terdekat"?

        Ini adalah inti dari seleksi state. Sistem melakukan langkah berikut secara berurutan:

        **Langkah 1 — Ambil pola target pada marker teramati**

        Dari marker yang sudah diketahui, kita bentuk sebuah pola:

        ```text
        Pola target (hanya marker teramati):
        {_obs_header}
        {_obs_values}
        ```

        **Langkah 2 — Bandingkan pola tiap reference terhadap pola target**

        Setiap reference haplotype diambil nilai allele-nya pada posisi yang sama,
        lalu dihitung berapa banyak yang cocok:

        ```text
        {_example_block}
        ```

        Semakin banyak posisi yang cocok, semakin "dekat" reference tersebut dengan target.

        **Langkah 3 — Urutkan reference dari yang paling mirip**

        Reference dengan pola paling cocok diurutkan ke atas. Kemudian diambil
        sejumlah reference teratas (diatur oleh slider **"Jumlah neighbour"** di sidebar)
        sebagai kandidat copying state untuk HMM.

        **Langkah 4 — Hasilnya: daftar state terpilih**

        Hanya reference yang terpilih inilah yang akan dipakai HMM. Reference lain
        yang polanya sangat berbeda dengan target tidak akan masuk hitungan.
        """
    )
    st.write("Tabel berikut menunjukkan detail kemiripan setiap reference terhadap target (diurutkan dari paling mirip):")
    st.dataframe(similarity_dataframe(reference_panel, target), width="stretch", hide_index=True)
    st.caption(
        "Kolom 'Jumlah cocok' menunjukkan berapa banyak posisi marker teramati yang nilainya sama "
        "antara reference dan target. Reference dengan jumlah cocok terbanyak cenderung masuk sebagai state HMM."
    )

with st.expander("Kalau PBWT asli di IMPUTE5, kira-kira seperti apa?"):
    st.markdown(
        """
        Demo ini memilih reference yang mirip dengan cara yang mudah dilihat: hitung berapa
        banyak allele yang cocok pada marker target yang sudah diketahui. PBWT asli lebih
        canggih dari itu.

        **Ide PBWT asli secara sederhana:**

        1. Pada setiap marker, reference haplotype diurutkan berdasarkan pola allele
           sebelumnya, dibaca dari marker saat ini ke arah kiri. Ini disebut **reverse prefix**.
        2. Haplotype yang berdekatan dalam urutan PBWT biasanya punya segmen allele yang
           panjang dan mirip.
        3. Target haplotype seolah-olah "disisipkan" ke dalam urutan itu.
        4. Reference haplotype yang berada di sekitar posisi target diambil sebagai kandidat state.

        Contoh kecil:

        ```text
        target pada marker teramati: M0=0, M1=1, M3=1, M5=0

        Jika sedang berada di M5, pola yang dilihat PBWT adalah dari kanan ke kiri:
        M5, M3, M1, M0 = 0, 1, 1, 0
        ```

        Reference yang memiliki pola kanan-ke-kiri mirip dengan `0,1,1,0` akan berada
        dekat dengan target dalam urutan PBWT. Karena itu reference tersebut dipilih
        sebagai copying state.

        **Perbedaan dengan demo ini:**

        | Bagian | Demo Streamlit | IMPUTE5 asli |
        |---|---|---|
        | Cara mencari kemiripan | hitung kecocokan sederhana | PBWT/FM-index dan neighbour/divergence selection |
        | Skala data | beberapa haplotype | ribuan sampai jutaan haplotype |
        | Tujuan | mudah dipahami | sangat cepat dan hemat memori |
        | Hasil seleksi | subset kecil reference | subset copying states untuk HMM |
        """
    )

st.subheader("3. HMM Yang Diperkecil")
st.write(
    "Versi sederhananya: HMM memberi bobot pada setiap reference terpilih. Reference yang lebih cocok "
    "dengan target akan mendapat bobot lebih besar. Bobot ini disebut posterior probability."
)
st.markdown(
    """
    Setelah seleksi state, ukuran HMM menjadi lebih kecil. HMM standar akan
    mempertimbangkan semua reference haplotype pada setiap marker. Di demo ini,
    hanya reference haplotype yang terpilih yang dipakai sebagai hidden state.

    - **Emission probability** mengecek apakah allele target cocok dengan allele reference.
    - **Transition probability** mengatur seberapa mudah model berpindah dari satu reference haplotype ke haplotype lain.
    - **Forward-backward** menggabungkan informasi dari sisi kiri dan kanan setiap marker.
    """
)
st.dataframe(target_display_dataframe(target), width="stretch")

with st.expander("Dari mana angka emission, transition, dan posterior muncul?", expanded=True):
    first_observed = next((index for index, allele in enumerate(target) if allele is not None), 0)

    st.markdown("### Ide besar HMM")
    st.markdown(
        'Bayangkan target haplotype seperti seseorang yang "menyalin" dari salah satu reference haplotype, '
        "satu marker demi satu marker. HMM bertugas menjawab:"
    )
    st.info("Pada marker ini, reference mana yang paling mungkin sedang disalin oleh target?")
    st.markdown(
        'Setiap reference terpilih disebut **hidden state**. Nama "hidden" karena kita tidak '
        "tahu pasti state mana yang benar — kita hanya bisa mengestimasi berdasarkan data."
    )

    st.divider()
    st.markdown("### 1. Emission probability — seberapa cocok allele-nya?")
    st.markdown(
        "Emission menjawab pertanyaan: *jika target sedang menyalin dari reference h_x, "
        "seberapa masuk akal allele target yang terlihat di marker ini?*"
    )
    st.markdown(
        "Aturannya sederhana:\n\n"
        "```text\n"
        "Jika allele target == allele reference → probabilitas = 1 - error_rate  (tinggi, cocok)\n"
        "Jika allele target != allele reference → probabilitas = error_rate       (rendah, tidak cocok)\n"
        "Jika target masih '?'                 → probabilitas = 1.0              (tidak ada info, netral)\n"
        "```"
    )
    st.markdown(f"Dengan **emission error rate** sekarang = `{error_rate:.3f}`:")
    st.markdown(
        f"```text\n"
        f"Cocok       → emission = {1.0 - error_rate:.3f}\n"
        f"Tidak cocok → emission = {error_rate:.3f}\n"
        f"```"
    )
    emission_lines = "\n".join(
        f"  h{ref_idx} allele di M{first_observed} = {reference_panel[ref_idx][first_observed]}"
        f" → {'cocok' if reference_panel[ref_idx][first_observed] == target[first_observed] else 'tidak cocok'}"
        f" → emission = {(1.0 - error_rate) if reference_panel[ref_idx][first_observed] == target[first_observed] else error_rate:.3f}"
        for ref_idx in selected_states
    )
    st.markdown(f"Contoh nyata dari data saat ini, pada marker **M{first_observed}** (nilai target = `{target[first_observed]}`):")
    st.markdown(f"```text\n{emission_lines}\n```")

    st.divider()
    st.markdown("### 2. Transition probability — seberapa mudah berpindah state?")
    st.markdown(
        "Transition menjawab: *apakah target tetap menyalin dari reference yang sama, "
        "atau berpindah ke reference lain di marker berikutnya?*\n\n"
        "Perpindahan ini merepresentasikan **rekombinasi** — fenomena biologis di mana "
        "dua kromosom bertukar segmen. Semakin tinggi switch rate, semakin sering "
        "model mempertimbangkan kemungkinan berpindah ke reference lain."
    )
    n_states = len(selected_states)
    denom = n_states - 1 if n_states > 1 else 1
    st.markdown(
        f"Dengan **recombination / switch rate** sekarang = `{recombination_rate:.2f}` "
        f"dan jumlah state terpilih = `{n_states}`:"
    )
    st.markdown(
        f"```text\n"
        f"P(tetap di state yang sama)  = 1 - switch_rate = {1.0 - recombination_rate:.2f}\n"
        f"P(pindah ke state lain)      = switch_rate / (jumlah state - 1)\n"
        f"                             = {recombination_rate:.2f} / {denom}\n"
        f"                             = {recombination_rate / denom:.3f}  (per state tujuan)\n"
        f"```"
    )
    st.markdown("Artinya, model cenderung \"setia\" pada satu reference, tapi tidak menutup kemungkinan pindah.")

    st.divider()
    st.markdown("### 3. Forward-backward — menggabungkan informasi kiri dan kanan")
    st.markdown(
        "Forward-backward adalah cara HMM menghitung **posterior probability** setiap state "
        "di setiap marker. Idenya:\n\n"
        "- **Forward pass** (kiri → kanan): di setiap marker, hitung seberapa besar probabilitas "
        "berada di state ini jika kita melihat semua allele dari awal hingga marker ini.\n"
        "- **Backward pass** (kanan → kiri): hitung hal yang sama tapi dari arah sebaliknya, "
        "yaitu semua allele dari akhir hingga marker ini.\n"
        "- **Posterior** = forward × backward, lalu dinormalisasi agar totalnya = 1.\n\n"
        "Hasilnya: setiap state mendapat \"bobot\" yang mencerminkan seberapa besar kontribusinya "
        "dalam menjelaskan pola allele target, dari **kedua sisi** sekaligus."
    )

    st.divider()
    st.markdown(f"Tabel berikut menunjukkan contoh angka emission dan transition untuk marker **M{first_observed}**:")
    st.dataframe(
        probability_example_dataframe(
            reference_panel,
            target,
            selected_states,
            first_observed,
            error_rate,
            recombination_rate,
        ),
        width="stretch",
        hide_index=True,
    )
    st.caption(
        f"Kolom 'Emission' menunjukkan seberapa cocok allele reference dengan target di M{first_observed}. "
        "Kolom 'P(tetap)' dan 'P(pindah)' menunjukkan peluang transisi ke marker berikutnya."
    )

with st.expander("Kalau HMM asli di IMPUTE5, kira-kira seperti apa?"):
    st.markdown(
        """
        IMPUTE5 memakai model probabilistik keluarga **Li-and-Stephens HMM**. Ide besarnya
        sama dengan demo ini: target haplotype dianggap menyalin potongan-potongan dari
        reference haplotype.

        **Hidden state asli**

        ```text
        state = reference haplotype yang sedang disalin oleh target
        ```

        Jika ada 10.000 reference haplotype, HMM standar bisa punya 10.000 state.
        IMPUTE5 mengurangi state ini memakai PBWT sebelum HMM dijalankan.

        **Emission probability asli**

        Emission menjawab pertanyaan:

        ```text
        Jika target sedang menyalin dari h3, seberapa masuk akal allele target yang terlihat?
        ```

        Jika allele target sama dengan allele reference, probabilitasnya tinggi.
        Jika berbeda, probabilitasnya rendah, tetapi tidak nol, karena bisa ada error,
        mutasi, atau ketidakcocokan model.

        **Transition probability asli**

        Transition menjawab pertanyaan:

        ```text
        Apakah target tetap menyalin dari reference yang sama,
        atau berpindah ke reference lain?
        ```

        Pada data genomik asli, peluang berpindah tidak hanya angka tetap seperti demo.
        Peluang ini dipengaruhi oleh **genetic map** dan jarak antar marker:

        - marker yang sangat dekat -> kemungkinan tetap di state yang sama lebih tinggi
        - marker yang lebih jauh -> kemungkinan berpindah state lebih besar

        **Forward-backward asli**

        Forward-backward menghitung probabilitas copying state di sepanjang kromosom.
        Setelah posterior state didapat, marker yang tidak ada pada target tetapi ada
        pada reference dapat diisi dengan cara:

        ```text
        P(allele = 1) = jumlah posterior state yang reference allele-nya 1
        ```

        Jadi angka probabilitas pada hasil imputasi berasal dari kombinasi:

        ```text
        kecocokan allele target-reference
        + peluang berpindah antar haplotype
        + informasi dari marker kiri dan kanan
        + allele reference pada marker yang hilang
        ```
        """
    )

with st.expander("Tampilkan posterior copying probability"):
    st.write(
        "Setiap baris adalah marker. Setiap kolom adalah reference haplotype yang terpilih. "
        "Nilai yang lebih besar berarti model lebih yakin target sedang menyalin dari haplotype tersebut."
    )
    st.dataframe(posterior_dataframe(selected_states, posterior), width="stretch")

st.subheader("4. Hasil Imputasi")
results = result_dataframe(target, imputed)
if results.empty:
    st.warning("Tidak ada marker `?` pada target, jadi tidak ada yang perlu diimputasi.")
else:
    st.write(
        "`P(allele = 1)` adalah keyakinan model bahwa marker yang hilang seharusnya bernilai `1`. "
        "Jika probabilitas ini di bawah 0.5, prediksi akhirnya menjadi `0`."
    )
    st.dataframe(results, width="stretch", hide_index=True)
    st.bar_chart(results.set_index("Marker")["P(allele = 1)"], width="stretch")

    explanation = explain_first_missing_marker(
        reference_panel,
        target,
        selected_states,
        posterior,
        imputed,
    )
    if explanation is not None:
        marker, explanation_table, probability = explanation
        with st.expander(f"Contoh cara membaca angka pada marker M{marker}", expanded=True):
            st.markdown(
                f"""
                Misalnya kita lihat marker `M{marker}`, salah satu marker yang awalnya `?`.

                Aplikasi melihat semua reference haplotype yang terpilih. Kalau sebuah reference
                punya allele `1` di `M{marker}`, maka bobot/posterior reference itu ikut dijumlahkan
                untuk menghitung `P(allele = 1)`.

                Kalau reference tersebut punya allele `0`, bobotnya tidak masuk ke penjumlahan
                `P(allele = 1)`.
                """
            )
            st.dataframe(explanation_table, width="stretch", hide_index=True)
            st.markdown(
                f"""
                Jadi untuk marker `M{marker}`:

                ```text
                P(allele = 1) = jumlah bobot reference terpilih yang allele-nya 1
                              = {probability:.3f}
                ```

                Jika nilainya lebih dari `0.5`, prediksi akhirnya menjadi `1`.
                Jika kurang dari `0.5`, prediksi akhirnya menjadi `0`.
                """
            )

    st.divider()
    st.subheader("Haplotype Target — Sebelum dan Sesudah Imputasi")
    st.write(
        "Di bawah ini adalah perbandingan target haplotype sebelum diimputasi (masih ada `?`) "
        "dan sesudah diimputasi (semua `?` sudah diganti dengan prediksi model)."
    )

    completed_values = []
    for i, allele in enumerate(target):
        if allele is not None:
            completed_values.append(str(allele))
        else:
            pred = 1 if imputed[i] >= 0.5 else 0
            completed_values.append(str(pred))

    before_row = ["?" if allele is None else str(allele) for allele in target]
    after_row = completed_values

    comparison_df = pd.DataFrame(
        [before_row, after_row],
        columns=marker_columns(len(target)),
        index=["Sebelum (dengan ?)", "Sesudah (imputasi)"],
    )
    st.dataframe(comparison_df, width="stretch")

    st.markdown("**Keterangan setiap marker:**")
    summary_rows = []
    for i, allele in enumerate(target):
        if allele is None:
            prob = imputed[i]
            pred = 1 if prob >= 0.5 else 0
            summary_rows.append({
                "Marker": f"M{i}",
                "Nilai asal": "?",
                "P(allele = 1)": f"{prob:.3f}",
                "Prediksi": pred,
                "Status": "diimputasi",
            })
        else:
            summary_rows.append({
                "Marker": f"M{i}",
                "Nilai asal": str(allele),
                "P(allele = 1)": "-",
                "Prediksi": allele,
                "Status": "sudah diketahui",
            })
    st.dataframe(pd.DataFrame(summary_rows), width="stretch", hide_index=True)
    st.caption(
        "Marker berstatus 'sudah diketahui' tidak diubah oleh model. "
        "Marker berstatus 'diimputasi' adalah tebakan model berdasarkan HMM."
    )

with st.expander("Penjelasan singkat alur implementasi"):
    st.markdown(
        """
        Implementasi ini mengikuti ide besar IMPUTE5:

        1. Mulai dari target haplotype yang memiliki marker hilang.
        2. Bandingkan target dengan reference panel yang lengkap.
        3. Pilih hanya reference haplotype yang mirip secara lokal sebelum komputasi HMM.
        4. Jalankan HMM yang lebih kecil dengan haplotype terpilih sebagai hidden state.
        5. Ubah posterior copying probability menjadi probabilitas allele untuk marker yang hilang.

        Poin pentingnya: IMPUTE5 tidak menghapus HMM. IMPUTE5 membuat HMM lebih
        scalable dengan mengurangi jumlah state sebelum forward-backward dijalankan.
        """
    )

with st.expander("Apa yang disederhanakan dibanding IMPUTE5 asli?"):
    st.markdown(
        """
        - Demo ini memakai array kecil berisi `0`, `1`, dan `?`, bukan file genomik nyata seperti VCF/BGEN/imp5.
        - Langkah PBWT di sini hanya PBWT-inspired; IMPUTE5 asli memakai PBWT/FM-index yang dioptimalkan.
        - Probabilitas HMM disederhanakan agar mudah dipahami.
        - IMPUTE5 asli menangani data phased diploid, genomic window, genetic map, dan reference panel besar.
        - IMPUTE5 asli dioptimalkan untuk kecepatan, memori, dan skala jutaan haplotype.
        """
    )

with st.expander("Bagaimana kalau memakai IMPUTE5 asli?"):
    st.markdown(
        """
        IMPUTE5 asli bukan library Python yang bisa langsung dipanggil seperti fungsi di demo ini.
        IMPUTE5 adalah software command-line bioinformatika, umumnya dijalankan di Linux atau WSL.

        Secara konsep, pipeline nyata terlihat seperti ini:

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
        impute5 \\
          --h reference_panel.bcf \\
          --g target_panel.bcf \\
          --m genetic_map.txt \\
          --r chr20:1000000-2000000 \\
          --o imputed_output.bgen
        ```

        Yang dibutuhkan jika memakai IMPUTE5 asli:

        - binary/software IMPUTE5,
        - data target yang biasanya sudah di-phase,
        - reference panel genomik,
        - genetic map,
        - format file yang sesuai seperti VCF/BCF/BGEN/imp5,
        - lingkungan Linux/WSL atau server/HPC.

        Untuk tugas presentasi, demo Streamlit ini lebih cocok untuk menjelaskan
        **bagaimana algoritmanya bekerja**. IMPUTE5 asli lebih cocok jika tujuannya
        adalah menjalankan pipeline bioinformatika nyata pada data genomik besar.
        """
    )