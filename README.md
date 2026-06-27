# LipinskyFilter - ADMETrics

A bioinformatics web application for analyzing molecular drug-likeness using **Lipinski's Rule of Five** and additional ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity) properties.

## 🧬 About

LipinskyFilter is a Streamlit-based tool designed for drug discovery researchers and chemists to quickly evaluate whether candidate compounds meet essential criteria for drug development. It calculates molecular descriptors including:

- **Lipinski's Rule of Five** - molecular weight, LogP, H-bond donors/acceptors
- **ADMET Properties** - solubility, bioavailability predictions
- **Molecular Visualization** - 2D and 3D structure rendering
- **Batch Analysis** - Process multiple compounds from CSV files

## ✨ Features

- **Interactive Molecular Input**: Enter SMILES notation or chemical structures
- **Real-time Property Calculation**: Instant calculation of drug-likeness metrics
- **3D Visualization**: Visualize molecular structures in 3D
- **Batch Processing**: Upload CSV files to analyze multiple compounds
- **Visual Dashboard**: Beautiful Plotly charts and metrics display
- **Dark Theme UI**: Modern, eye-friendly interface with cyan accents

## 🛠️ Tech Stack

- **[Streamlit](https://streamlit.io/)** - Web framework for data apps
- **[RDKit](https://www.rdkit.org/)** - Cheminformatics library
- **[Plotly](https://plotly.com/)** - Interactive visualizations
- **[Pandas](https://pandas.pydata.org/)** - Data manipulation
- **[py3Dmol](https://3dmol.csb.pitt.edu/)** - 3D molecular visualization
- **[OpenCV](https://opencv.org/)** - Image processing

## 🚀 Getting Started

### Prerequisites

- Python 3.8+
- pip or conda

### Installation

1. **Clone the repository**
   ```bash
   git clone <repository-url>
   cd LipinskyFilter
   ```

2. **Create a virtual environment** (optional but recommended)
   ```bash
   python -m venv venv
   # On Windows:
   venv\Scripts\activate
   # On macOS/Linux:
   source venv/bin/activate
   ```

3. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

### Running Locally

```bash
streamlit run app.py
```

The application will start at `http://localhost:8501`

## 📊 Usage

### Single Compound Analysis
1. Enter a SMILES string (e.g., `CC(C)Cc1ccc(cc1)C(C)C(=O)O` for Ibuprofen)
2. View calculated properties and Lipinski compliance
3. Visualize the 2D/3D structure

### Batch Analysis
1. Upload a CSV file with a SMILES column
2. The tool processes all compounds
3. Download results as CSV with all calculated properties

### CSV Format
```
SMILES,Name
CC(C)Cc1ccc(cc1)C(C)C(=O)O,Ibuprofen
CC(=O)Oc1ccccc1C(=O)O,Aspirin
```

## 📋 Project Structure

```
LipinskyFilter/
├── app.py              # Main Streamlit application
├── compounds.csv       # Sample compound data
├── requirements.txt    # Python dependencies
├── vercel.json         # Vercel deployment configuration
├── venv/              # Virtual environment
└── README.md          # This file
```

## 🌐 Deployment

### Deploy to Vercel

The project is configured for easy deployment to Vercel:

```bash
vercel
```

The `vercel.json` configuration handles routing and Python runtime setup.

### Environment Variables

No environment variables required for basic functionality.

## 📚 Lipinski's Rule of Five

A compound meets Lipinski's criteria if it has:
- **Molecular Weight** ≤ 500 Da
- **LogP** ≤ 5
- **H-bond Donors** ≤ 5
- **H-bond Acceptors** ≤ 10

Drugs satisfying these criteria typically have better oral bioavailability.

## 🔧 Development

### Adding New Features

1. Modify `app.py` to add new properties or visualizations
2. Update `requirements.txt` if adding new dependencies
3. Test locally before deployment

### Common Modifications

- **Change theme colors**: Edit CSS in `st.markdown()` section
- **Add new molecular descriptors**: Use RDKit's `Descriptors` module
- **Modify page layout**: Adjust Streamlit columns and containers

## 🐛 Troubleshooting

**Issue**: RDKit installation fails
- **Solution**: Use conda instead: `conda install -c conda-forge rdkit`

**Issue**: Molecules not displaying
- **Solution**: Ensure SMILES strings are valid using RDKit validation

**Issue**: Streamlit not finding modules
- **Solution**: Verify virtual environment is activated

## 📝 License

[Specify your license here, e.g., MIT, GPL, Apache 2.0]

## 👨‍💻 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -m 'Add improvement'`)
4. Push to branch (`git push origin feature/improvement`)
5. Open a Pull Request

## 📧 Support

For issues, questions, or feature requests, please open an issue on the GitHub repository.

## 🙏 Acknowledgments

- RDKit community for cheminformatics tools
- Streamlit for the web framework
- Plotly for interactive visualizations
