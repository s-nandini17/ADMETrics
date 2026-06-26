import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import os
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, AllChem
from stmol import showmol
import py3Dmol

# ==========================================
# 1. PAGE CONFIGURATION
# ==========================================
st.set_page_config(
    page_title="ADMETrics",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed"
)

# ==========================================
# 2. SYNTAX THEME CSS (Slate, Cyan, Purple)
# ==========================================
st.markdown("""
<style>
    /* Dark Slate Backgrounds */
    .main {
        background-color: #0F111A;
        color: #E2E8F0;
    }
    
    h1, h2, h3, h4, p, span, div {
        color: #E2E8F0;
        font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    }

    /* Metric Cards */
    .metric-card {
        background-color: #1E1E2E;
        padding: 24px;
        border-radius: 8px;
        border: 1px solid #333344;
        border-top: 4px solid #00E5FF; /* Cyan Accent */
        box-shadow: 0 4px 6px rgba(0,0,0,0.3);
        text-align: center;
        margin-bottom: 1rem;
    }
    
    .metric-value {
        font-size: 28px;
        font-weight: 700;
        color: #FFFFFF;
        margin-bottom: 4px;
    }
    
    .metric-label {
        font-size: 13px;
        text-transform: uppercase;
        letter-spacing: 0.05em;
        color: #94A3B8;
    }

    /* Inspector Violation Tags */
    .violation-tag {
        display: inline-block;
        background-color: rgba(217, 70, 239, 0.15); /* Purple Tint */
        color: #D946EF; /* Electric Purple */
        border: 1px solid #D946EF;
        padding: 6px 12px;
        border-radius: 6px;
        margin: 4px;
        font-weight: 600;
        font-size: 14px;
    }
    .pass-tag {
        display: inline-block;
        background-color: rgba(0, 229, 255, 0.15); /* Cyan Tint */
        color: #00E5FF; 
        border: 1px solid #00E5FF;
        padding: 6px 12px;
        border-radius: 6px;
        margin: 4px;
        font-weight: 600;
        font-size: 14px;
    }

    /* Hide default Streamlit styling */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    
    /* Input and Selectbox styling */
    .stSelectbox > div > div {
        background-color: #1E1E2E;
        color: white;
    }
</style>
""", unsafe_allow_html=True)

# ==========================================
# 3. CORE LOGIC FUNCTIONS
# ==========================================
@st.cache_data
def calculate_lipinski_properties(smiles):
    """Calculate properties from SMILES"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return None
        return {
            'MW': round(Descriptors.MolWt(mol), 2),
            'LogP': round(Crippen.MolLogP(mol), 2),
            'HBD': int(Descriptors.NumHDonors(mol)),
            'HBA': int(Descriptors.NumHAcceptors(mol))
        }
    except:
        return None

def get_violation_details(props):
    """Return explicit list of violated rules"""
    if props is None: return []
    violations = []
    if props['MW'] > 500: violations.append(f"MW > 500 ({props['MW']} Da)")
    if props['LogP'] > 5: violations.append(f"LogP > 5 ({props['LogP']})")
    if props['HBD'] > 5: violations.append(f"HBD > 5 ({props['HBD']})")
    if props['HBA'] > 10: violations.append(f"HBA > 10 ({props['HBA']})")
    return violations

def get_status_label(violation_count):
    if violation_count == 0: return "Pass"
    elif violation_count == 1: return "1 Violation"
    else: return "2+ Violations"

def generate_3d_mol(smiles):
    """Convert SMILES to 3D MolBlock for py3Dmol"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
        return Chem.MolToMolBlock(mol)
    except:
        return None

# ==========================================
# 4. APP LAYOUT & ARCHITECTURE
# ==========================================
st.markdown("<h1>🧬 ADMETrics</h1>", unsafe_allow_html=True)
st.markdown("<p style='color: #94A3B8; margin-bottom: 2rem;'>High-Throughput Molecular Property & Drug-Likeness Dashboard</p>", unsafe_allow_html=True)

# Define the 3 Pages (Tabs)
tab1, tab2, tab3 = st.tabs(["📊 Overview & Upload", "📈 System Visualizations", "🔬 Inspector & Export"])

# --- TAB 1: OVERVIEW & UPLOAD ---
with tab1:
    st.markdown("""
    ### 📋 Lipinski Criteria
    Evaluating oral bioavailability based on the Rule of Five:
    * **Molecular Weight (MW)** ≤ 500 Da
    * **Lipophilicity (LogP)** ≤ 5
    * **H-bond Donors (HBD)** ≤ 5
    * **H-bond Acceptors (HBA)** ≤ 10
    ---
    """)
    
    col_upload, col_default = st.columns([3, 1])
    with col_upload:
        uploaded_file = st.file_uploader("Upload CSV (columns: compound_id, name, smiles)", type="csv")
    with col_default:
        st.markdown("<br>", unsafe_allow_html=True) 
        default_csv = "compounds.csv"
        if os.path.exists(default_csv) and st.button("📂 Load Default Dataset"):
            uploaded_file = default_csv

# ==========================================
# 5. DATA PROCESSING & VISUALIZATION
# ==========================================
if uploaded_file is not None or os.path.exists("compounds.csv"):
    try:
        df = pd.read_csv(uploaded_file)
        
        with st.spinner("Processing molecular structures..."):
            df['Properties'] = df['smiles'].apply(calculate_lipinski_properties)
            df_valid = df[df['Properties'].notna()].copy()
            
            if len(df_valid) == 0:
                st.error("No valid SMILES strings found.")
                st.stop()
                
            # Expand properties
            df_valid['MW'] = df_valid['Properties'].apply(lambda x: x['MW'])
            df_valid['LogP'] = df_valid['Properties'].apply(lambda x: x['LogP'])
            df_valid['HBD'] = df_valid['Properties'].apply(lambda x: x['HBD'])
            df_valid['HBA'] = df_valid['Properties'].apply(lambda x: x['HBA'])
            df_valid['Violation_Details'] = df_valid['Properties'].apply(get_violation_details)
            df_valid['Violation_Count'] = df_valid['Violation_Details'].apply(len)
            df_valid['Status'] = df_valid['Violation_Count'].apply(get_status_label)
            
            # Format display names for dropdowns
            df_valid['display_name'] = df_valid['compound_id'].astype(str) + " - " + df_valid['name']

        # --- TAB 1 (Continued): METRICS ---
        with tab1:
            st.markdown("### Batch Summary")
            pass_count = (df_valid['Violation_Count'] == 0).sum()
            one_viol = (df_valid['Violation_Count'] == 1).sum()
            multi_viol = (df_valid['Violation_Count'] >= 2).sum()
            
            m1, m2, m3, m4 = st.columns(4)
            m1.markdown(f'<div class="metric-card"><div class="metric-value">{len(df_valid)}</div><div class="metric-label">Total Valid</div></div>', unsafe_allow_html=True)
            m2.markdown(f'<div class="metric-card" style="border-top-color: #00E5FF;"><div class="metric-value" style="color:#00E5FF;">{pass_count}</div><div class="metric-label">Pass (0)</div></div>', unsafe_allow_html=True)
            m3.markdown(f'<div class="metric-card" style="border-top-color: #A78BFA;"><div class="metric-value" style="color:#A78BFA;">{one_viol}</div><div class="metric-label">Caution (1)</div></div>', unsafe_allow_html=True)
            m4.markdown(f'<div class="metric-card" style="border-top-color: #D946EF;"><div class="metric-value" style="color:#D946EF;">{multi_viol}</div><div class="metric-label">Fail (≥2)</div></div>', unsafe_allow_html=True)

        # --- TAB 2: SYSTEM VISUALIZATIONS ---
        with tab2:
            st.markdown("### Molecular Weight vs LogP")
            
            # Syntax Theme Colors mapping
            color_map = {'Pass': '#00E5FF', '1 Violation': '#A78BFA', '2+ Violations': '#D946EF'}
            
            fig = px.scatter(
                df_valid, x='MW', y='LogP', color='Status', 
                hover_data=['name', 'HBD', 'HBA'],
                color_discrete_map=color_map,
                template='plotly_dark',
            )
            
            # Threshold lines
            fig.add_hline(y=5, line_dash="dash", line_color="#475569")
            fig.add_vline(x=500, line_dash="dash", line_color="#475569")
            
            fig.update_layout(
                height=450, margin=dict(t=20),
                paper_bgcolor='#0F111A', plot_bgcolor='#1E1E2E'
            )
            st.plotly_chart(fig, use_container_width=True)
            
            st.markdown("### Global Distributions")
            c1, c2 = st.columns(2)
            with c1:
                fig_mw = px.histogram(df_valid, x='MW', nbins=30, template='plotly_dark', color_discrete_sequence=['#00E5FF'])
                fig_mw.add_vline(x=500, line_dash="solid", line_color="#D946EF")
                fig_mw.update_layout(height=300, margin=dict(t=20), paper_bgcolor='#0F111A', plot_bgcolor='#1E1E2E')
                st.plotly_chart(fig_mw, use_container_width=True)
            with c2:
                fig_logp = px.histogram(df_valid, x='LogP', nbins=30, template='plotly_dark', color_discrete_sequence=['#A78BFA'])
                fig_logp.add_vline(x=5, line_dash="solid", line_color="#D946EF")
                fig_logp.update_layout(height=300, margin=dict(t=20), paper_bgcolor='#0F111A', plot_bgcolor='#1E1E2E')
                st.plotly_chart(fig_logp, use_container_width=True)

        # --- TAB 3: INSPECTOR & EXPORT ---
        with tab3:
            st.markdown("### 🔬 Compound Inspector")
            
            # Dropdown to select a compound
            selected_comp = st.selectbox("Select a compound for deep-dive analysis:", df_valid['display_name'].tolist())
            
            if selected_comp:
                # Isolate selected data
                target_data = df_valid[df_valid['display_name'] == selected_comp].iloc[0]
                
                col_3d, col_radar = st.columns([1, 1])
                
                with col_3d:
                    st.markdown(f"**3D Structure:** {target_data['name']}")
                    mol_block = generate_3d_mol(target_data['smiles'])
                    if mol_block:
                        view = py3Dmol.view(width=400, height=350)
                        view.addModel(mol_block, 'mol')
                        view.setStyle({'stick': {'colorscheme': 'cyanCarbon'}})
                        view.setBackgroundColor('#1E1E2E')
                        view.zoomTo()
                        showmol(view, height=350, width=400)
                    else:
                        st.warning("Could not generate 3D coordinates for this SMILES.")

                with col_radar:
                    st.markdown("**Rule of Five Adherence** (Scaled to Limits)")
                    
                    # Scale values relative to Lipinski thresholds (Threshold = 1.0)
                    radar_vals = [
                        target_data['MW'] / 500,
                        max(target_data['LogP'], 0) / 5, # Avoid negative logP breaking the chart
                        target_data['HBD'] / 5,
                        target_data['HBA'] / 10
                    ]
                    categories = ['MW (Lim: 500)', 'LogP (Lim: 5)', 'HBD (Lim: 5)', 'HBA (Lim: 10)']
                    
                    fig_radar = go.Figure()
                    
                    # Add Threshold Boundary (Perfect Hexagon at 1.0)
                    fig_radar.add_trace(go.Scatterpolar(
                        r=[1, 1, 1, 1], theta=categories, fill='toself',
                        name='Lipinski Threshold', line_color='rgba(255, 255, 255, 0.3)',
                        fillcolor='rgba(255, 255, 255, 0.05)'
                    ))
                    
                    # Add Compound Data
                    line_color = '#00E5FF' if target_data['Violation_Count'] == 0 else '#D946EF'
                    fig_radar.add_trace(go.Scatterpolar(
                        r=radar_vals, theta=categories, fill='toself',
                        name=target_data['name'], line_color=line_color,
                        fillcolor=line_color.replace(')', ', 0.4)').replace('rgb', 'rgba')
                    ))
                    
                    fig_radar.update_layout(
                        polar=dict(radialaxis=dict(visible=True, range=[0, max(1.5, max(radar_vals))])),
                        showlegend=False, height=280, margin=dict(t=20, b=20, l=40, r=40),
                        paper_bgcolor='#0F111A'
                    )
                    st.plotly_chart(fig_radar, use_container_width=True)
                    
                    # Explicit Violation Display
                    st.markdown("**Status & Violations:**")
                    if target_data['Violation_Count'] == 0:
                        st.markdown("<span class='pass-tag'>✓ Fully Compliant</span>", unsafe_allow_html=True)
                    else:
                        for v in target_data['Violation_Details']:
                            st.markdown(f"<span class='violation-tag'>⚠️ {v}</span>", unsafe_allow_html=True)

            st.markdown("---")
            st.markdown("### 💾 Dataset Explorer & Export")
            
            # Simple clean dataframe for export view
            display_cols = ['compound_id', 'name', 'MW', 'LogP', 'HBD', 'HBA', 'Violation_Count', 'Status']
            st.dataframe(df_valid[display_cols], use_container_width=True, height=250, hide_index=True)
            
            st.download_button(
                label="📥 Export Dataset to CSV",
                data=df_valid[display_cols].to_csv(index=False),
                file_name="admetrics_results.csv",
                mime="text/csv",
                type="primary"
            )

    except Exception as e:
        st.error(f"Error processing file: {str(e)}")