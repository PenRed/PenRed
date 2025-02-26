
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
        f.write(f"materials/{name}/eabs/gamma {gammaCutoffValue}\n")
    else:
        f.write(f"materials/{name}/range/gamma {gammaCutoffValue}\n")

    # Electron
    if electronCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/electron {electronCutoffValue}\n")
    else:
        f.write(f"materials/{name}/range/electron {electronCutoffValue}\n")

    # Positron
    if positronCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/positron {positronCutoffValue}\n")
    else:
        f.write(f"materials/{name}/range/positron {positronCutoffValue}\n")

    if advanced:
        f.write(f"materials/{name}/C1 {C1}\n")
        f.write(f"materials/{name}/C2 {C2}\n")
        f.write(f"materials/{name}/WCC {WCC}\n")
        f.write(f"materials/{name}/WCR {WCR}\n")

        
    f.write(f"materials/{name}/density {density}\n")
    for e in composition:
        f.write(f"materials/{name}/elements/{e[0]} {e[1]}\n")
    f.write("\n")
