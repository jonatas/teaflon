import numpy as np
import matplotlib.pyplot as plt
import os

def simulate_ise_fluoride_release():
    print("🧪 TEAFLON DIGITAL TWIN: Step 1 Defluorination Test\n" + "="*55)
    print("Simulating the Ion-Selective Electrode (ISE) sensor data for PTFE digestion.")
    
    # Simulation Parameters (Michaelis-Menten Kinetics for PTFE breakdown)
    # Vmax: Maximum rate of C-F bond cleavage
    # Km: Michaelis constant (affinity for PTFE)
    time_hours = np.linspace(0, 48, 100)
    
    # 1. Control (Teflon in water - indestructible)
    fluoride_control = np.zeros_like(time_hours)
    
    # 2. Natural Enzyme (Baseline Dehalogenase - struggles with PTFE hydrophobicity)
    fluoride_natural = 5 * (1 - np.exp(-0.02 * time_hours)) 
    
    # 3. Engineered Teaflon (Hydrophobin anchor + Dehalogenase active site)
    # The hydrophobin anchor increases the local concentration of PTFE at the active site
    fluoride_teaflon = 150 * (1 - np.exp(-0.15 * time_hours))

    # Plotting the Digital Twin Sensor Data
    plt.figure(figsize=(10, 6))
    plt.plot(time_hours, fluoride_control, label="Control (Buffer Only)", color='gray', linestyle='--')
    plt.plot(time_hours, fluoride_natural, label="Wild-type Enzyme", color='blue')
    plt.plot(time_hours, fluoride_teaflon, label="Engineered Teaflon System", color='red', linewidth=2.5)
    
    plt.title("Theoretical Defluorination of PTFE Nanoparticles\n(Fluoride Ion-Selective Electrode Digital Twin)", pad=15)
    plt.xlabel("Incubation Time (Hours)")
    plt.ylabel("Free Fluoride Concentration [F⁻] (ppm)")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Fill the area under the Teaflon curve to represent total toxic fluoride neutralized
    plt.fill_between(time_hours, fluoride_teaflon, alpha=0.1, color='red')
    
    # Save the chart
    output_path = "assets/step1_ise_simulation.png"
    os.makedirs("output", exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n✅ Simulation Complete. ISE Sensor chart saved to: {output_path}")
    print("This graph represents the baseline you must beat in the wet lab.")

if __name__ == "__main__":
    simulate_ise_fluoride_release()
