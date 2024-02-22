# %%
import uuid
import random
import pandas as pd
from q2_moshpit.busco.utils import _draw_busco_plots_for_render
import altair as alt

# Enable all rows
alt.data_transformers.disable_max_rows()


# Function to generate random fractions
def generate_fractions(n):
    fractions = sorted(random.random() for _ in range(n - 1))
    fractions = [fractions[0]] + [
        fractions[i] - fractions[i - 1] for i in range(1, n - 1)
        ] + [1 - fractions[-1]]
    return fractions


# Generate data
data = []
for i in range(1, 2):  # samples, brakes at 66 * 20 (still ok)
    for j in range(1, 66*24):  # mags. with unlimit 66 * 24 is still ok
        unique_id = uuid.uuid4()
        numbers = generate_fractions(5)
        data.append(
            [
                f"{unique_id}.fasta",
                "bacteria_odb10",
                numbers[0],
                numbers[1],
                numbers[2],
                numbers[3],
                numbers[4],
                124,
                170295,
                170295,
                "0.000%",
                27,
                f"sample{i}"
            ]
        )

# Define column names
columns = [
    "Input_file",
    "Dataset",
    "Complete",
    "Single",
    "Duplicated",
    "Fragmented",
    "Missing",
    "n_markers",
    "Scaffold N50",
    "Contigs N50",
    "Percent gaps",
    "Number of scaffolds",
    "sample_id"
]

# Create DataFrame
df = pd.DataFrame(data, columns=columns)

# Draw
dictto = _draw_busco_plots_for_render(
    df,
    width=600,
    height=30,
    titleFontSize=20,
    labelFontSize=17,
    spacing=20
)

dictto
# %%
