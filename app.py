import streamlit as st
from Bio import SeqIO
from CRISPR import CRISPR
import matplotlib.pyplot as plt
import random
from typing import Dict, Any

# ============================================================
# DATA & UTILITY FUNCTIONS
# ============================================================

# Define file paths
VARIANTS_FILES = {
    "Wuhan (Original COVID strain)": "wuhan.fasta",
    "Alpha (UK Variant)": "alpha.fasta",
    "Delta (Indian Variant)": "delta.fasta",
    "Omicron (South African Variant)": "omicron.fasta"
}

@st.cache_data(show_spinner="Loading Spike Gene Sequences...")
def load_sequence(filename: str):
    """Loads and caches a single spike gene sequence from a FASTA file."""
    try:
        record = SeqIO.read(filename, "fasta")
        return str(record.seq)
    except FileNotFoundError:
        st.error(f"Error: FASTA file not found ({filename}).")
        return None
    except Exception as e:
        st.error(f"Error reading {filename}: {e}")
        return None

@st.cache_data(show_spinner="Loading all sequences for comparison...")
def load_all_sequences(variant_files: Dict[str, str]) -> Dict[str, str]:
    """Loads and caches all spike gene sequences for the comparison tab."""
    sequences = {}
    for name_long, filename in variant_files.items():
        name_short = name_long.split('(')[0].strip()
        sequences[name_short] = load_sequence(filename)
    return {k: v for k, v in sequences.items() if v is not None}


def create_outcome_pie_chart(result: Dict[str, Any]):
    """Generates the Matplotlib pie chart for simulation outcomes."""
    off_risk = float(result.get("total_offtarget_risk", 0))
    p_cleave = float(result.get("p_cleave", 0))
    
    # Simplified probabilities for educational display
    loss_func = p_cleave * 0.4
    partial_func = p_cleave * 0.3
    unmodeled = 0.01
    no_effect = max(0, 1 - (off_risk + loss_func + partial_func + unmodeled))
    
    labels = [
        "Off-Target Cuts",
        "Loss of Function",
        "Partial Function",
        "No Significant Effect",
        "Unmodeled"
    ]
    values = [off_risk, loss_func, partial_func, no_effect, unmodeled]
    colors = ["#4B8BBE", "#FF6F61", "#FFD166", "#06D6A0", "#BDBDBD"]

    fig, ax = plt.subplots(figsize=(5, 5))
    wedges, texts, autotexts = ax.pie(
        values,
        startangle=140,
        colors=colors,
        autopct=lambda p: f'{p:.1f}%' if p >= 2.0 else '',
        textprops={'fontsize': 9, 'color': 'black'}
    )

    ax.legend(
        wedges,
        labels,
        title="Outcome Types",
        loc="center left",
        bbox_to_anchor=(1, 0, 0.5, 1),
        fontsize=9
    )

    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(9)
        autotext.set_fontweight('bold')

    ax.set_title("Possible Outcomes of CRISPR Edit", fontsize=12, fontweight='bold', pad=15)
    ax.axis("equal")
    plt.tight_layout()
    return fig


# ============================================================
# TAB 1 LOGIC: CRISPR SIMULATION
# ============================================================

def run_crispr_simulation_tab(variant_files: Dict[str, str]):
    """Handles the UI and logic for the CRISPR Simulation tab."""
    
    with st.sidebar.form("crispr_settings_form"):
        st.header("🔬 Experiment Settings")

        # Variant Choice
        variant_choice_key = st.selectbox(
            "🦠 Choose which COVID variant to explore",
            list(variant_files.keys()),
            key="variant_choice_sim"
        )
        
        genome_seq = load_sequence(variant_files[variant_choice_key])
        if not genome_seq:
            return 
        
        st.write(f"📏 Spike gene length: **{len(genome_seq)} base pairs**")

        # Position Slider (1-based for user)
        max_pos = len(genome_seq) - 21
        if max_pos < 1:
            max_pos = 1
        default_position = min(1000, max_pos)

        position = st.slider(
            "🔍 Where should CRISPR cut the gene? (1-based position shown)",
            1, max_pos, default_position,
            help="Slide to choose where CRISPR should make its cut."
        )

        # Convert to 0-based for internal indexing
        pos0 = position - 1

        # Safe guide extraction (20 bases)
        guide_start = min(pos0, max(0, len(genome_seq) - 20))
        guide_end = min(len(genome_seq), guide_start + 20)
        guide_seq = genome_seq[guide_start:guide_end]

        if len(guide_seq) < 20:
            st.warning("Selected position is near the end; guide is shorter than 20 bases.")

        st.markdown(f"**🧬 CRISPR Guide (20 bases):** `{guide_seq}`")

        mutation_type = st.radio(
            "🧪 What type of change do you want to make?",
            ["Tiny Deletion", "Small Insertion", "Tiny Mutation (swap letters)"],
            help="Choose the type of edit CRISPR will make (conceptual)."
        )

        submitted = st.form_submit_button("🚀 Run CRISPR Simulation")

    if submitted:
        st.subheader(f"🧫 Simulation Results for {variant_choice_key}")

        try:
            crispr_model = CRISPR()
            result = crispr_model.run(guide_seq, genome_seq)
        except NameError:
            st.error("CRISPR class not found. Ensure the 'CRISPR' module is available.")
            return
        except Exception as e:
            st.error(f"CRISPR model failed to run: {e}")
            return

        repair_labels = {
            "frameshift": "Big Change (Loss of Function)",
            "small_del": "Small Deletion",
            "ins_1bp": "Tiny Insertion",
            "inframe": "Small Change (Function Maintained)",
            "unknown": "No Major Change"
        }
        repair_text = repair_labels.get(result.get("repair_outcome"), "Minor DNA Edit")

        # Display metrics
        p_cleave = float(result.get("p_cleave", 0))
        off_risk = float(result.get("total_offtarget_risk", 0))

        col1, col2, col3 = st.columns(3)
        col1.metric("✂️ Cut Probability", f"{p_cleave:.2f}")
        col2.metric("🧬 Repair Outcome", repair_text)
        col3.metric("⚠️ Off-Target Risk", f"{off_risk:.2f}")

        # Off-target table
        if result.get("offtargets") and len(result["offtargets"]) > 0:
            st.markdown("### ⚠️ Possible Off-Target Sites (Top 5)")
            st.dataframe(result["offtargets"][:5])
        else:
            st.markdown("✅ No strong off-target sites predicted.")

        # Pie chart
        st.markdown("### 📊 Predicted Outcome Distribution")
        fig = create_outcome_pie_chart(result)
        st.pyplot(fig)

        # Gene snippet visualization
        start = max(0, pos0 - 10)
        end = min(len(genome_seq), pos0 + 30)
        snippet = genome_seq[start:end]
        guide_start_idx = pos0 - start
        if guide_start_idx < 0:
            guide_start_idx = 0

        st.markdown("### 🧫 Zoom into the Spike Gene Around the Cut Site")
        st.code(f"Genome Snippet:\n{snippet}\n{' ' * guide_start_idx}{'^' * 5} ← CRISPR Cut")

        st.markdown("""
        ---
        ### 🧠 What this means:
        - **Cut Probability:** Likelihood that Cas9 cuts at the selected spot.  
        - **Off-Target Risk:** Chance CRISPR cuts elsewhere in the genome.  
        - **Repair Outcome:** What kind of mutation happens after repair.  
        - **Pie Chart:** Estimated frequency of each possible edit type.
        """)

    else:
        st.info("""
👈 Use the sidebar to:
1. Pick a COVID variant  
2. Move the slider to select a cut site  
3. Pick the type of edit  
4. Then click **🚀 Run CRISPR Simulation**!
""")


# ============================================================
# TAB 2 LOGIC: COMPARE VARIANTS
# ============================================================

def run_compare_variants_tab(variant_files: Dict[str, str]):
    """Handles the UI and logic for the Compare Variants tab."""
    
    sequences = load_all_sequences(variant_files)
    if not sequences:
        return 

    st.subheader("🧬 Compare Spike Gene Sequences Across Variants")
    st.markdown("""
    This section highlights how the **Spike gene** differs between variants.  
    Each colored letter shows a change from the **Wuhan (original)** sequence.
    """)

    base_key = "Wuhan"
    base = sequences.get(base_key, "")
    
    compare_options = [k for k in sequences.keys() if k != base_key]
    if not compare_options:
        st.warning("Sequence data is not available for comparison.")
        return

    compare_choice = st.selectbox("🧬 Choose a variant to compare with Wuhan:", compare_options)
    target = sequences.get(compare_choice, "")

    if not base or not target:
        st.warning(f"Sequence data missing for {base_key} or {compare_choice}.")
        return

    compare_length = min(500, len(base), len(target))
    highlighted = []
    for i in range(compare_length):
        if base[i] != target[i]:
            highlighted.append(f"<span style='color:red; font-weight:bold;'>{target[i]}</span>")
        else:
            highlighted.append(target[i])
    highlighted_html = "".join(highlighted)
    if len(base) > compare_length:
        highlighted_html += "..."

    st.markdown(f"**Comparing {base_key} vs {compare_choice} (First {compare_length} bases):**")
    st.markdown(
        f"<pre style='white-space: pre-wrap; word-break: break-all; font-size: 14px;'>{highlighted_html}</pre>",
        unsafe_allow_html=True
    )
    st.caption("Red, bold letters = mutation compared to Wuhan sequence.")

    st.markdown("""
    ---
    💡 **Try this:** Select different variants above to see how their Spike gene changed —  
    that’s what made Omicron so different and harder for vaccines to target.
    """)


# ============================================================
# MAIN APP FLOW
# ============================================================

st.set_page_config(
    page_title="🧬 CRISPR Spike Gene Simulator",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 CRISPR Spike Gene Simulator")

try:
    st.image("crispr.png", caption="CRISPR–Cas9 complex binding to DNA")
except Exception as e:
    st.warning(f"⚠️ Could not load the image: {e}")



st.markdown("""
Explore how **CRISPR gene editing** could work on the *Spike gene* of different COVID-19 variants.  
This interactive demo makes complex biology simple for students and educators.
""")

tab1, tab2 = st.tabs(["🧫 CRISPR Simulation", "🧬 Compare Variants"])
with tab1:
    run_crispr_simulation_tab(VARIANTS_FILES)
with tab2:
    run_compare_variants_tab(VARIANTS_FILES)

st.markdown("""
---
🔬 *Educational demo for exploring CRISPR and SARS-CoV-2 Spike gene edits.*
""")
