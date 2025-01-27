# Master Thesis Project: Computational Neuroscience

This repository contains the code for my **master's thesis project** in computational neuroscience. The project involves an **interactive simulation** of the learning process in a **cerebellum module**, focusing on its main constituents: granule cells, olive cells, and a Purkinje cell.

## Theoretical Background

The theoretical foundation of this project is best described in the following articles:

1. **Shakirov, V., Altunina, O., Shaposhnikov, D., et al.**  
   *Cerebellar plasticity-based equalization of total input to inferior olive cells.*  
   *Neurosci Behav Physi*, 53, 739–751 (2023).  
   [https://doi.org/10.1007/s11055-023-01422-8](https://doi.org/10.1007/s11055-023-01422-8)

2. **Shakirov, V., Dorofeev, V., & Dunin-Barkowski, W.**  
   *Cerebellar plasticity-based equalization of total input to inferior olive cells: properties of the model dynamics.*  
   *Neurosci Behav Physi*, 53, 729–738 (2023).  
   [https://doi.org/10.1007/s11055-023-01423-7](https://doi.org/10.1007/s11055-023-01423-7)

3. **Shakirov, V., Dorofeev, V., Lebedev, A., et al.**  
   *Dynamic chaos in cerebellum and electrical synapses between climbing fiber cells of inferior olives.*  
   *Neurosci Behav Physi*, 53, 717–728 (2023).  
   [https://doi.org/10.1007/s11055-023-01420-w](https://doi.org/10.1007/s11055-023-01420-w)

## Project Details

- **Development Environment**: This project was created using **Microsoft Visual C++ (MFC)**.  
- **Main Working File**: The core functionality of the project can be found in `ProbeDlg.cpp`.  

The code simulates how cerebellar plasticity processes affect learning dynamics in the context of input equalization for inferior olive cells. This implementation aims to provide insights into the role of specific neural components and their interactions within the cerebellum.

## Code Usage

### Prerequisites
Before building and running the project, ensure the following prerequisites are met:
- **Development Environment**:  
  Microsoft Visual Studio (recommended: **Visual Studio 2008**).  
- **Target Platform**:  
  Windows 32-bit (**Win32**) architecture.
- **Libraries and Dependencies**:  
  - MFC (Microsoft Foundation Classes) must be installed and enabled in Visual Studio. 

### Compilation Instructions
1. **Open the Project**:
   - Clone the repository or download the project files.
   - Open the `probe.sln` file in **Visual Studio**.

2. **Configure Build Settings**:
   - Ensure the configuration is set to `Debug|Win32` or `Release|Win32` (as defined in the project file).

3. **Build the Project**:
   - From the Visual Studio menu, select **Build → Build Solution** or press `Ctrl+Shift+B`.
   - The executable (`probe.exe`) will be created in the `.\Debug` or `.\Release` directory.

### Notes
- **Language Encoding**:  
  The project files use `windows-1251` encoding, suitable for Cyrillic characters. Ensure this encoding is preserved when editing the files.
- **Compatibility**:  
  The project uses older Visual Studio formats (`.dsp`, `.dsw`, and `.vcproj`). For newer Visual Studio versions, the project may require conversion or compatibility adjustments.

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
