
def createMaterial(f, name, index,
                   gammaCutoffType, gammaCutoffValue,
                   electronCutoffType, electronCutoffValue,
                   positronCutoffType, positronCutoffValue,
                   density, composition,
                   advanced,
                   C1, C2, WCC, WCR):
    f.write(f"# Material '{name}' configuration\n")    
    f.write(f"materials/{name}/number {index}\n")
    f.write(f"materials/{name}/filename \"{name}.mat\"\n")
    
    # Gamma
    if gammaCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/gamma {gammaCutoffValue:.3e}\n")
    else:
        f.write(f"materials/{name}/range/gamma {gammaCutoffValue:.3e}\n")

    # Electron
    if electronCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/electron {electronCutoffValue:.3e}\n")
    else:
        f.write(f"materials/{name}/range/electron {electronCutoffValue:.3e}\n")

    # Positron
    if positronCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/positron {positronCutoffValue:.3e}\n")
    else:
        f.write(f"materials/{name}/range/positron {positronCutoffValue:.3e}\n")

    if advanced:
        f.write(f"materials/{name}/C1 {C1:.2f}\n")
        f.write(f"materials/{name}/C2 {C2:.2f}\n")
        f.write(f"materials/{name}/WCC {WCC:.3e}\n")
        f.write(f"materials/{name}/WCR {WCR:.3e}\n")

        
    f.write(f"materials/{name}/density {density:.4e}\n")
    for e in composition:
        f.write(f"materials/{name}/elements/{e[0]} {e[1]:.4e}\n")
    f.write("\n")
