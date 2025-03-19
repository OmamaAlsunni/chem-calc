import sys 


class Element:
  def __init__(self,symbol,name,atomic_number,atomic_mass,group,period):
    self.name=name
    self.symbol=symbol
    self.atomic_number=atomic_number
    self.atomic_mass=atomic_mass
    self.group=group
    self.period=period
  def describe(self):
    return f"Element: {self.name} (Symbol: {self.symbol}) has an atomic number of {self.atomic_number}, an atomic mass of {self.atomic_mass}, is in group {self.group}, and is in period {self.period}"


# Main Group Elements
hydrogen = Element("H", "Hydrogen", 1, 1.008, "1", "1")
helium = Element("He", "Helium", 2, 4.003, "18", "1")
lithium = Element("Li", "Lithium", 3, 6.941, "1", "2")
beryllium = Element("Be", "Beryllium", 4, 9.012, "2", "2")
boron = Element("B", "Boron", 5, 10.811, "13", "2")
carbon = Element("C", "Carbon", 6, 12.011, "14", "2")
nitrogen = Element("N", "Nitrogen", 7, 14.007, "15", "2")
oxygen = Element("O", "Oxygen", 8, 15.999, "16", "2")
fluorine = Element("F", "Fluorine", 9, 18.998, "17", "2")
neon = Element("Ne", "Neon", 10, 20.180, "18", "2")
sodium = Element("Na", "Sodium", 11, 22.990, "1", "3")
magnesium = Element("Mg", "Magnesium", 12, 24.305, "2", "3")
aluminum = Element("Al", "Aluminum", 13, 26.982, "13", "3")
silicon = Element("Si", "Silicon", 14, 28.085, "14", "3")
phosphorus = Element("P", "Phosphorus", 15, 30.974, "15", "3")
sulfur = Element("S", "Sulfur", 16, 32.065, "16", "3")
chlorine = Element("Cl", "Chlorine", 17, 35.453, "17", "3")
argon = Element("Ar", "Argon", 18, 39.948, "18", "3")
potassium = Element("K", "Potassium", 19, 39.098, "1", "4")
calcium = Element("Ca", "Calcium", 20, 40.078, "2", "4")

# Transition Elements
scandium = Element("Sc", "Scandium", 21, 44.956, "3", "4")
titanium = Element("Ti", "Titanium", 22, 47.867, "4", "4")
vanadium = Element("V", "Vanadium", 23, 50.942, "5", "4")
chromium = Element("Cr", "Chromium", 24, 51.996, "6", "4")
manganese = Element("Mn", "Manganese", 25, 54.938, "7", "4")
iron = Element("Fe", "Iron", 26, 55.845, "8", "4")
cobalt = Element("Co", "Cobalt", 27, 58.933, "9", "4")
nickel = Element("Ni", "Nickel", 28, 58.693, "10", "4")
copper = Element("Cu", "Copper", 29, 63.546, "11", "4")
zinc = Element("Zn", "Zinc", 30, 65.380, "12", "4")
gallium = Element("Ga", "Gallium", 31, 69.723, "13", "4")
germanium = Element("Ge", "Germanium", 32, 72.640, "14", "4")
arsenic = Element("As", "Arsenic", 33, 74.922, "15", "4")
selenium = Element("Se", "Selenium", 34, 78.960, "16", "4")
bromine = Element("Br", "Bromine", 35, 79.904, "17", "4")
krypton = Element("Kr", "Krypton", 36, 83.798, "18", "4")

rubidium = Element("Rb","Rubidium",37,85.468,"1","5")
strontium = Element("Sr","Strontium",38,87.62,"2","5")
cesium = Element("Cs","Cesium",55,132.91,"1","6")
barium = Element("Ba","Barium",56,137.33,"2","6")
yttrium = Element("Y","Yttrium",39,88.906,"3","5")
zirconium = Element("Zr","Zirconium",40,91.224,"4","5")
niobium = Element("Nb","Niobium",41,92.906,"5","5")
molybdenum = Element("Mo","Molybdenum",42,95.95,"6","5")
technetium = Element("Tc","Technetium",43,98,"7","5")
ruthenium = Element("Ru","Ruthenium",44,101.07,"8","5")
rhodium = Element("Rh","Rhodium",45,102.91,"9","5")
palladium = Element("Pd","Palladium",46,106.42,"10","5")
silver = Element("Ag","Silver",47,107.8682,"11","5")
cadmium = Element("Cd","Cadmium",48,112.41,"12","5")
indium = Element("In","Indium",49,114.82,"13","5")
tin = Element("Sn","Tin",50,118.71,"14","5")
antimony = Element("Sb","Antimony",51,121.76,"15","5")
tellurium = Element("Te","Tellurium",52,127.6,"16","5")
iodine = Element("I","Iodine",53,126.90,"17","5")
xenon = Element("Xe","Xenon",54,131.29,"18","5")
lanthanum = Element("La","Lanthanum",57,138.91,"3","6")
cerium = Element("Ce","Cerium",58,140.12,"4","6")
praseodymium = Element("Pr","Praseodymium",59,140.91,"5","6")
neodymium = Element("Nd","Neodymium",60,144.24,"6","6")
promethium = Element("Pm","Promethium",61,145,"7","6")
samarium = Element("Sm","Samarium",62,150.36,"7","6")
europium = Element("Eu","Europium",63,151.96,"7","6")
gadolinium = Element("Gd","Gadolinium",64,157.25,"7","6")
terbium = Element("Tb","Terbium",65,158.93,"7","6")
dysprosium = Element("Dy","Dysprosium",66,162.50,"7","6")
holmium = Element("Ho","Holmium",67,164.93,"7","6")
erbium = Element("Er","Erbium",68,167.26,"7","6")
thulium = Element("Tm","Thulium",69,168.93,"7","6")
ytterbium = Element("Yb","Ytterbium",70,173.05,"7","6")
lutetium = Element("Lu","Lutetium",71,174.97,"7","6")
hafnium = Element("Hf","Hafnium",72,178.49,"4","6")
tantalum = Element("Ta","Tantalum",73,180.95,"5","6")
tungsten = Element("W","Tungsten",74,183.84,"6","6")
rhenium = Element("Re","Rhenium",75,186.21,"7","6")
osmium = Element("Os","Osmium",76,190.23,"8","6")
iridium = Element("Ir","Iridium",77,192.22,"9","6")
platinum = Element("Pt","Platinum",78,195.08,"10","6")
gold = Element("Au","Gold",79,196.9665,"11","6")
mercury = Element("Hg","Mercury",80,200.59,"12","6")
thallium = Element("Tl","Thallium",81,204.38,"13","6")
lead = Element("Pb","Lead",82,207.2,"14","6")
bismuth = Element("Bi","Bismuth",83,208.98,"15","6")
polonium = Element("Po","Polonium",84,209,"16","6")
astatine = Element("At","Astatine",85,210,"17","6")
radon = Element("Rn","Radon",86,222,"18","6")
francium = Element("Fr","Francium",87,223,"1","7")
radium = Element("Ra","Radium",88,226,"2","7")
actinium = Element("Ac","Actinium",89,227,"3","7")
thorium = Element("Th","Thorium",90,232.04,"4","7")
protactinium = Element("Pa","Protactinium",91,231.04,"5","7")
uranium = Element("U","Uranium",92,238.03,"6","7")
neptunium = Element("Np","Neptunium",93,237,"7","7")
plutonium = Element("Pu","Plutonium",94,244,"7","7")
americium = Element("Am","Americium",95,243,"7","7")
curium = Element("Cm","Curium",96,247,"7","7")
berkelium = Element("Bk","Berkelium",97,247,"7","7")
californium = Element("Cf","Californium",98,251,"7","7")
einsteinium = Element("Es","Einsteinium",99,252,"7","7")
fermium = Element("Fm","Fermium",100,257,"7","7")
mendelevium = Element("Md","Mendelevium",101,258,"7","7")
nobelium = Element("No","Nobelium",102,259,"7","7")
lawrencium = Element("Lr","Lawrencium",103,262,"7","7")
rutherfordium = Element("Rf","Rutherfordium",104,267,"7","7")
dubnium = Element("Db","Dubnium",105,268,"7","7")
seaborgium = Element("Sg","Seaborgium",106,269,"7","7")
bohrium = Element("Bh","Bohrium",107,270,"7","7")
hassium = Element("Hs","Hassium",108,277,"7","7")
meitnerium = Element("Mt","Meitnerium",109,278,"7","7")
darmstadtium = Element("Ds","Darmstadtium",110,281,"7","7")
roentgenium = Element("Rg","Roentgenium",111,282,"7","7")
copernicium = Element("Cn","Copernicium",112,285,"7","7")
nihonium = Element("Nh","Nihonium",113,286,"7","7")
flerovium = Element("Fl","Flerovium",114,289,"7","7")
moscovium = Element("Mc","Moscovium",115,290,"7","7")
livermorium = Element("Lv","Livermorium",116,293,"7","7")
tennessine = Element("Ts","Tennessine",117,294,"7","7")
oganesson = Element("Og","Oganesson",118,294,"7","7")




elements = {
    # Main Group Elements
    "H": hydrogen, "Hydrogen": hydrogen,
    "He": helium, "Helium": helium,
    "Li": lithium, "Lithium": lithium,
    "Be": beryllium, "Beryllium": beryllium,
    "B": boron, "Boron": boron,
    "C": carbon, "Carbon": carbon,
    "N": nitrogen, "Nitrogen": nitrogen,
    "O": oxygen, "Oxygen": oxygen,
    "F": fluorine, "Fluorine": fluorine,
    "Ne": neon, "Neon": neon,
    "Na": sodium, "Sodium": sodium,
    "Mg": magnesium, "Magnesium": magnesium,
    "Al": aluminum, "Aluminum": aluminum,
    "Si": silicon, "Silicon": silicon,
    "P": phosphorus, "Phosphorus": phosphorus,
    "S": sulfur, "Sulfur": sulfur,
    "Cl": chlorine, "Chlorine": chlorine,
    "Ar": argon, "Argon": argon,
    "K": potassium, "Potassium": potassium,
    "Ca": calcium, "Calcium": calcium,

    # Transition Elements
    "Sc": scandium, "Scandium": scandium,
    "Ti": titanium, "Titanium": titanium,
    "V": vanadium, "Vanadium": vanadium,
    "Cr": chromium, "Chromium": chromium,
    "Mn": manganese, "Manganese": manganese,
    "Fe": iron, "Iron": iron,
    "Co": cobalt, "Cobalt": cobalt,
    "Ni": nickel, "Nickel": nickel,
    "Cu": copper, "Copper": copper,
    "Zn": zinc, "Zinc": zinc,
    "Y": yttrium, "Yttrium": yttrium,
    "Zr": zirconium, "Zirconium": zirconium,
    "Nb": niobium, "Niobium": niobium,
    "Mo": molybdenum, "Molybdenum": molybdenum,
    "Tc": technetium, "Technetium": technetium,
    "Ru": ruthenium, "Ruthenium": ruthenium,
    "Rh": rhodium, "Rhodium": rhodium,
    "Pd": palladium, "Palladium": palladium,
    "Ag": silver, "Silver": silver,
    "Cd": cadmium, "Cadmium": cadmium,
    "Hf": hafnium, "Hafnium": hafnium,
    "Ta": tantalum, "Tantalum": tantalum,
    "W": tungsten, "Tungsten": tungsten,
    "Re": rhenium, "Rhenium": rhenium,
    "Os": osmium, "Osmium": osmium,
    "Ir": iridium, "Iridium": iridium,
    "Pt": platinum, "Platinum": platinum,
    "Au": gold, "Gold": gold,
    "Hg": mercury, "Mercury": mercury,

    # Post-Transition Elements
    "Ga": gallium, "Gallium": gallium,
    "Ge": germanium, "Germanium": germanium,
    "As": arsenic, "Arsenic": arsenic,
    "Se": selenium, "Selenium": selenium,
    "Br": bromine, "Bromine": bromine,
    "Kr": krypton, "Krypton": krypton,
    "In": indium, "Indium": indium,
    "Sn": tin, "Tin": tin,
    "Sb": antimony, "Antimony": antimony,
    "Te": tellurium, "Tellurium": tellurium,
    "I": iodine, "Iodine": iodine,
    "Xe": xenon, "Xenon": xenon,
    "Tl": thallium, "Thallium": thallium,
    "Pb": lead, "Lead": lead,
    "Bi": bismuth, "Bismuth": bismuth,
    "Po": polonium, "Polonium": polonium,
    "At": astatine, "Astatine": astatine,
    "Rn": radon, "Radon": radon,

    # Alkali and Alkaline Earth Metals
    "Rb": rubidium, "Rubidium": rubidium,
    "Sr": strontium, "Strontium": strontium,
    "Cs": cesium, "Cesium": cesium,
    "Ba": barium, "Barium": barium,
    "Fr": francium, "Francium": francium,
    "Ra": radium, "Radium": radium,

    # Lanthanides
    "La": lanthanum, "Lanthanum": lanthanum,
    "Ce": cerium, "Cerium": cerium,
    "Pr": praseodymium, "Praseodymium": praseodymium,
    "Nd": neodymium, "Neodymium": neodymium,
    "Pm": promethium, "Promethium": promethium,
    "Sm": samarium, "Samarium": samarium,
    "Eu": europium, "Europium": europium,
    "Gd": gadolinium, "Gadolinium": gadolinium,
    "Tb": terbium, "Terbium": terbium,
    "Dy": dysprosium, "Dysprosium": dysprosium,
    "Ho": holmium, "Holmium": holmium,
    "Er": erbium, "Erbium": erbium,
    "Tm": thulium, "Thulium": thulium,
    "Yb": ytterbium, "Ytterbium": ytterbium,
    "Lu": lutetium, "Lutetium": lutetium,

    # Actinides
    "Ac": actinium, "Actinium": actinium,
    "Th": thorium, "Thorium": thorium,
    "Pa": protactinium, "Protactinium": protactinium,
    "U": uranium, "Uranium": uranium,
    "Np": neptunium, "Neptunium": neptunium,
    "Pu": plutonium, "Plutonium": plutonium,
    "Am": americium, "Americium": americium,
    "Cm": curium, "Curium": curium,
    "Bk": berkelium, "Berkelium": berkelium,
    "Cf": californium, "Californium": californium,
    "Es": einsteinium, "Einsteinium": einsteinium,
    "Fm": fermium, "Fermium": fermium,
    "Md": mendelevium, "Mendelevium": mendelevium,
    "No": nobelium, "Nobelium": nobelium,
    "Lr": lawrencium, "Lawrencium": lawrencium,

    # Superheavy Elements
    "Rf": rutherfordium, "Rutherfordium": rutherfordium,
    "Db": dubnium, "Dubnium": dubnium,
    "Sg": seaborgium, "Seaborgium": seaborgium,
    "Bh": bohrium, "Bohrium": bohrium,
    "Hs": hassium, "Hassium": hassium,
    "Mt": meitnerium, "Meitnerium": meitnerium,
    "Ds": darmstadtium, "Darmstadtium": darmstadtium,
    "Rg": roentgenium, "Roentgenium": roentgenium,
    "Cn": copernicium, "Copernicium": copernicium,
    "Nh": nihonium, "Nihonium": nihonium,
    "Fl": flerovium, "Flerovium": flerovium,
    "Mc": moscovium, "Moscovium": moscovium,
    "Lv": livermorium, "Livermorium": livermorium,
    "Ts": tennessine, "Tennessine": tennessine,
    "Og": oganesson, "Oganesson": oganesson
}

from math import gcd

# ==================== 1. MOLE & MASS FUNCTIONS ====================

def mole_mass():
  while True:
    print("\nMole & mass calculations")
    print("1.find moles")
    print("2.find mass")
    print("3.find compound molar mass")
    print("4.element molar mass")
    print("5.return")
    mole_choice = input("\nenter your choice:").strip()
    if mole_choice == "1":
      find_moles()
    elif mole_choice == "2":
      find_mass()
    elif mole_choice == "3":
      compound_molar_mass()
    elif mole_choice == "4":
      find_element_molar_mass()
    elif mole_choice == "5":
      break
    else:
      print("\ninvalid entry")

def find_moles():
  print("\nFindng number of moles:")
  try:
    mass = input("mass(g):").strip()
    mass = float(mass)
    element = input("write the element's name or symbol..").strip()
    element = element.capitalize()
    if element in elements:
      element_obj = elements[element]
      moles = mass/element_obj.atomic_mass
      print(f"\nNumber of moles: {moles} mol")
      return moles
    else:
      print("\nElement not found. Please provide its molar mass.")
      molar_mass = input("Molar mass(g/mol):").strip()
      molar_mass = float(molar_mass)
      moles = mass/molar_mass
      print(f"\nNumber of moles: {moles} mol")
      return moles
  except ValueError:
    print("\nError: Please enter valid numeric values for mass and/or molar mass.")
    return None
  except ZeroDivisionError:
    print("\nError: Molar mass cannot be zero.")
    return None

def find_mass():
  print("\nFindng mass:")
  try:
    moles = input("Number of moles(mol):").strip()
    moles = float(moles)
    element = input("write the element's name or symbol..").strip()
    element = element.capitalize()
    if element in elements:
      element_obj = elements[element]
      mass = moles*element_obj.atomic_mass
      print(f"\nMass: {mass} g")
      return mass
    else:
      print("Element not found. Please provide its molar mass.")
      molar_mass = input("Molar mass(g/mol):").strip()
      molar_mass = float(molar_mass)
      mass = moles*molar_mass
      print(f"Mass: {mass} g")
      return mass
  except ValueError:
    print("\nError: Please enter valid numeric values.")
    return None

def find_element_molar_mass():
  element = input("search for an element you want to know its molar mass(use name or symbol):").strip().capitalize()
  if element in elements:
    element_obj = elements[element]
    print(f"\n{element} molar mass: {element_obj.atomic_mass} g/mol")
  else:
    print(f"{element} does not correspond to any currently recognized chemical element in the periodic table")

def compound_molar_mass():
  compound_molar_mass = 0
  more_elements = "Yes"
  while more_elements == "Yes":
    element = input("write the element's name or symbol:").strip().capitalize()
    if element in elements:
      element_obj = elements[element]
      try:
        atoms = int(input(f"number of atoms of {element} in the compound:").strip())
        compound_molar_mass += atoms*element_obj.atomic_mass
      except ValueError:
        print("Invalid input. Please enter a numeric value for the number of atoms.")
        continue
    else:
      print(f"{element} does not correspond to any currently recognized chemical element in the periodic table.")
      choice = input("continue by providing the molar mass of the element?(Yes/No):").strip().capitalize()
      if choice == "Yes":
        try:
          molar_mass = float(input("Molar mass(g/mol):").strip())
          atoms = int(input("number of atoms of the element in the compound:").strip())
          compound_molar_mass += atoms*molar_mass
        except ValueError:
          print("Invalid input. Please enter a numeric value for the molar mass and the number of atoms")
          continue
      elif choice == "No":
        continue
    more_elements = input("does the compound consist of any more elements?(Yes/No):").strip().capitalize()
  print(f"\nCompound molar mass: {compound_molar_mass} g/mol")
  return compound_molar_mass

# ==================== 2. STOICHIOMETRY FUNCTIONS ====================

def stoichiometry():
  while True:
    print("\nStoichiometry")
    print("1.Empirical & molecular formulas")
    print("2.Formulas of hydrates")
    print("3.return")
    stoichiometry_choice = input("\nenter your choice:").strip()
    if stoichiometry_choice == "1":
      empirical_formula_calculation()
    elif stoichiometry_choice == "2":
      calculate_hydrate_formula()
    elif stoichiometry_choice == "3":
      break
    else:
      print("\ninvalid entry")

def find_smallest_multiplier(values):
  decimals = [val%1 for val in values if val%1 != 0]
  if not decimals:
    return 1
  denominators = [round(1/d) for d in decimals]
  lcm = denominators[0]
  for num in denominators[1:]:
    lcm = (lcm*num)//gcd(lcm,num)
  return lcm

def empirical_formula_calculation():
  moles_dict = {}
  more_elements = "Yes"
  while more_elements == "Yes":
    element = input("write the element's name or symbol:").strip().capitalize()
    if element in elements:
      element_obj = elements[element]
      try:
        element_mass = float(input(f"mass of {element} in the sample (g):"))
        element_moles = element_mass/element_obj.atomic_mass
        moles_dict[element] = element_moles
      except ValueError:
        print(f"Invalid input. Please enter a valid number for the mass of {element}.")
        continue
    else:
      print(f"{element} does not correspond to any currently recognized chemical element in the periodic table")
      continue
    more_elements = input("does the sample consist of any more elements?(Yes/No):").strip().capitalize()

  min_moles = min(moles_dict.values())
  ratio_dict = {elem:mol/min_moles for elem,mol in moles_dict.items()}

  multiplier = find_smallest_multiplier(ratio_dict.values())
  final_ratios = {elem: round(num*multiplier) for elem,num in ratio_dict.items()}

  subscript_digits = str.maketrans("0123456789","₀₁₂₃₄₅₆₇₈₉")
  empirical_formula = "".join([f"{elem}{str(num).translate(subscript_digits) if num>1 else''}" for elem,num in final_ratios.items()])
  print(f"\nEmpirical formula: {empirical_formula}")

  return empirical_formula

def calculate_hydrate_formula(hydrated_salt_mass=None, anhydrous_salt_mass=None, anhydrous_salt_molar_mass=None):
  if hydrated_salt_mass is None and anhydrous_salt_mass is None:
    hydrated_salt_mass = float(input("Enter the mass of the hydrated salt (g): "))
    anhydrous_salt_mass = float(input("Enter the mass of the anhydrous salt (g): "))
  anhydrous_salt_molar_mass = compound_molar_mass()

  water_mass = hydrated_salt_mass - anhydrous_salt_mass
  n_water = water_mass/18.015
  n_salt = anhydrous_salt_mass/anhydrous_salt_molar_mass
  x = n_water/n_salt
  print(f"\nChemical formula: Salt ⋅ {round(x)}H₂O")

# ==================== 3. SOLUTION CONCENTRATION FUNCTIONS ====================

def solution_concentration():
  while True:
    print("\nSolution concentration")
    print("1.molarity\n2.molality\n3.mole fraction\n4.mass percentage\n5.volume percentage\n6.return")
    concentration_choice = input("\nenter your choice:")
    if concentration_choice == "1":
      molarity_calculations()
    elif concentration_choice == "2":
      molality_calculations()
    elif concentration_choice == "3":
      mole_fraction_calculation()
    elif concentration_choice == "4":
      mass_percentage_calculation()
    elif concentration_choice == "5":
      volume_percentage_calculation()
    elif concentration_choice == "6":
      break
    else:
      print("invalid entry")

def molarity_calculations():
  while True:
    print("\nMolarity calculations")
    print("1.Find molarity")
    print("2.Find moles")
    print("3.find volume")
    print("4.return")
    molarity_choice = input("\nenter your choice:")
    if molarity_choice == "1":
      find_molarity()
    elif molarity_choice == "2":
      find_moles_from_molarity()
    elif molarity_choice == "3":
      find_volume_from_molarity()
    elif molarity_choice == "4":
      break
    else:
      print("\ninvalid entry")

def find_molarity():
  try:
    moles = input("solute moles(mol):")
    volume = input("solution volume(L):")
    moles = float(moles)
    volume = float(volume)

    if volume == 0:
      print("\nError: Volume cannot be zero.")
      return None

    molarity = moles/volume
    print("\nmolarity:", molarity, "M")
    return molarity
  except ValueError:
    print("\nError: Please enter valid numeric values.")
    return None
  except ZeroDivisionError:
    print("\nError: Volume cannot be zero.")
    return None

def find_moles_from_molarity():
  molarity = input("\nmolarity:")
  volume = input("solution volume(L):")
  print("\nmoles:", float(molarity)*float(volume), "mol")

def find_volume_from_molarity():
  molarity = input("\nmolarity:")
  moles = input("solute moles(mol):")
  print("\nvolume:", float(moles)/float(molarity), "L")

def molality_calculations():
  while True:
    print("\nMolality calculations\n1.find molality\n2.find moles\n3.find mass\n4.return")
    molality_choice = input("\nenter your choice:")
    if molality_choice == "1":
      find_molality()
    elif molality_choice == "2":
      find_moles_from_molality()
    elif molality_choice == "3":
      find_mass_from_molality()
    elif molality_choice == "4":
      break
    else:
      print("invalid entry")

def find_molality():
  moles = input("solute moles(mol):")
  mass = input("solvent mass(kg):")
  molality = float(moles)/float(mass)
  print("molality:", molality, "mol/kg")

def find_moles_from_molality():
  molality = input("molality(mol/kg):")
  mass = input("solvent mass(kg):")
  moles = float(molality)*float(mass)
  print("solute moles:", moles, "mol")

def find_mass_from_molality():
  molality = input("molality(mol/kg):")
  moles = input("solute moles(mol):")
  mass = float(moles)/float(molality)
  print("solvent mass:", mass, "kg")

def mole_fraction_calculation():
  print("\nMole fraction")
  molesA = input("solute moles(mol):")
  molesB = input("solution moles(mol):")
  mole_fraction = float(molesA)/float(molesB)
  print("\nmole fraction:", mole_fraction)

def mass_percentage_calculation():
  print("\nMass percentage")
  massA = input("solute mass(kg):")
  massB = input("solution mass(kg):")
  mass_percentage = float(massA)*100/float(massB)
  print("\nmass percentage:", mass_percentage, "%")

def volume_percentage_calculation():
  print("\nVolume percentage")
  volumeA = input("solute volume(L):")
  volumeB = input("solution volume(L):")
  volume_percentage = float(volumeA)*100/float(volumeB)
  print("\nvolume percentage:", volume_percentage, "%")

# ==================== 4. UNIT CONVERSION FUNCTIONS ====================

def unit_conversions():
  while True:
    print("\nUnit conversions")
    print("1.moles - particles")
    print("2.liters - moles")
    print("3.temperature conversions")
    print("4.energy conversions")
    print("5.pressure conversions")
    print("6.return")
    unit_choice = input("\nenter your choice:")
    if unit_choice == "1":
      moles_particles_conversion()
    elif unit_choice == "2":
      liters_moles_conversion()
    elif unit_choice == "3":
      temperature_conversion()
    elif unit_choice == "4":
      energy_conversion()
    elif unit_choice == "5":
      pressure_conversion()
    elif unit_choice == "6":
      break
    else:
      print("\ninvalid entry") 

def moles_particles_conversion():
  while True:
    print("\nMoles - particles")
    print("1.find number of  particles")
    print("2.find number of moles")
    print("3.return")
    choice = input("\nenter your choice:")
    if choice == "1":
      find_particles_from_moles()
    elif choice == "2":
      find_moles_from_particles()
    elif choice == "3":
      break
    else:
      print("\ninvalid entry")

def find_particles_from_moles():
  try:
    moles = input("number of moles(mol):")
    moles = float(moles)
    avogadro = 6.0221367*10**23
    print("\nparticles:", moles*avogadro, "particle(s)")
  except ValueError:
    print("\nError: Please enter a valid number for moles.")

def find_moles_from_particles():
  particles = input("number of particles(e+power):").strip()
  avogadro = 6.0221367*10**23
  print("\nmoles:", float(particles)/float(avogadro), "mol")

def liters_moles_conversion():
  while True:
    print("\nLiters - moles")
    print("1.find number of moles")
    print("2.find number of liters")
    print("3.return")
    choice = input("\nenter your choice:")
    if choice == "1":
      find_moles_from_volume()
    elif choice == "2":
      find_volume_from_moles()
    elif choice == "3":
      break
    else:
      print("\ninvalid entry")

def find_moles_from_volume():
  volume = input("volume(L):")
  print("\nmoles:", float(volume)/22.4, "mol")

def find_volume_from_moles():
  moles = input("number of moles(mol):")
  print("\nliters:", float(moles)*22.4, "L")

def temperature_conversion():
  while True:
    print("\nTemperature conversions")
    print("1.kelvin to celsius")
    print("2.celsius to kelvin")
    print("3.return")
    choice = input("\nenter your choice:")
    if choice == "1":
      kelvin_to_celsius()
    elif choice == "2":
      celsius_to_kelvin()
    elif choice == "3":
      break
    else:
      print("\ninvalid entry")

def kelvin_to_celsius():
  kelvin = input("temperature in kelvin:")
  celsius = float(kelvin)-273.15
  print("\ntemperature in celsius:", celsius, "°C")

def celsius_to_kelvin():
  celsius = input("temperature in celsius:")
  kelvin = float(celsius)+273.15
  print("\ntemperature in kelvin:", kelvin, "K")

def energy_conversion():
  while True:
    print("\nEnergy conversions")
    print("1.joules to calories(Cal)")
    print("2.calories to joules(J)")
    print("3.return")
    choice = input("\nenter your choice:")
    if choice == "1":
      joules_to_calories()
    elif choice == "2":
      calories_to_joules()
    elif choice == "3":
      break
    else:
      print("\ninvalid entry")

def joules_to_calories():
  joules = input("energy in joules(J):")
  calories = float(joules)/4.184
  print("\nenergy in calories:", calories, "Cal")

def calories_to_joules():
  calories = input("energy in calories(Cal):")
  joules = float(calories)*4.184
  print("\nenergy in joules:", joules, "J")

def pressure_conversion():
  while True:
    print("\nPressure conversions")
    print("1.atm to Pa")
    print("2.Pa to atm")
    print("3.mmHg to Pa")
    print("4.Pa to mmHg")
    print("5.atm to mmHg")
    print("6.mmHg to atm")
    print("7.return")
    choice = input("\nenter your choice:")
    if choice == "1":
      atm_to_pa()
    elif choice == "2":
      pa_to_atm()
    elif choice == "3":
      mmhg_to_pa()
    elif choice == "4":
      pa_to_mmhg()
    elif choice == "5":
      atm_to_mmhg()
    elif choice == "6":
      mmhg_to_atm()
    elif choice == "7":
      break
    else:
      print("\ninvalid entry")

def atm_to_pa():
  atm = input("pressure in atm:")
  pa = float(atm)*101325
  print("\npreasure in Pa:", pa, "Pa")

def pa_to_atm():
  pa = input("pressure in Pa:")
  atm = float(pa)/101325
  print("\npressure in atm:", atm, "atm")

def mmhg_to_pa():
  mmHg = input("pressure in mmHg:")
  pa = float(mmHg)*133.322
  print("\npressure in Pa:", pa, "Pa")

def pa_to_mmhg():
  pa = input("pressure in Pa:")
  mmHg = float(pa)/133.322
  print("\npressure in mmHg:", mmHg, "mmHg")

def atm_to_mmhg():
  atm = input("pressure in atm:")
  mmHg = float(atm)*760
  print("\npressure in mmHg:", mmHg, "mmHg")

def mmhg_to_atm():
  mmHg = input("pressure in mmHg:")
  atm = float(mmHg)/760
  print("\npressure in atm:", atm, "atm")

# ==================== 5. PERCENT YIELD FUNCTION ====================

def percent_yield():
  try:
    actual_yield = input("enter the actual yield:")
    theoretical_yield = input("enter the theoretical yield:")
    actual_yield = float(actual_yield)
    theoretical_yield = float(theoretical_yield)

    if theoretical_yield == 0:
      print("\nError: Theoretical yield cannot be zero.")
      return None

    percentage = actual_yield/theoretical_yield*100
    print("\npercent yield:", percentage, "%")
    return percentage
  except ValueError:
    print("\nError: Please enter valid numeric values for yields.")
    return None
  except ZeroDivisionError:
    print("\nError: Theoretical yield cannot be zero.")
    return None

# ==================== 6. PERCENTAGE PER MASS FUNCTION ====================

def percentage_per_mass():
  results = []
  loop = "Yes"
  while loop == "Yes":
    print("\nPercentage per mass")
    try:
      element_mass = input("element mass:")
      compound_mass = input("compund mass:")
      element_mass = float(element_mass)
      compound_mass = float(compound_mass)

      if compound_mass == 0:
        print("\nError: Compound mass cannot be zero.")
        continue

      percentage = element_mass*100/compound_mass
      print("\npercentage per mass:", percentage, "%")
      results.append(percentage)
      loop = input("\nwould you like to calculate this again(Yes/No)?").strip().capitalize()
    except ValueError:
      print("\nError: Please enter valid numeric values.")
      loop = input("\nwould you like to try again(Yes/No)?").strip().capitalize()
    except ZeroDivisionError:
      print("\nError: Compound mass cannot be zero.")
      loop = input("\nwould you like to try again(Yes/No)?").strip().capitalize()

  return results

# ==================== EXIT FUNCTION ====================

def exit_calculator():
  print("\nThank you for using this calculator! Have a nice day!")
  sys.exit()

# ==================== MAIN MENU FUNCTION ====================

def show_menu():
  while True:
    print("\nchemical calculator")
    print("===================")
    print("1.mole & mass calculations")
    print("2.stoichiometry")
    print("3.solution concentration")
    print("4.unit conversions")
    print("5.percent yield calculation")
    print("6.percentage per mass calculation")
    print("7.exit")
    choice = input("\nenter your choice:")
    if choice == "1":
      mole_mass()
    elif choice == "2":
      stoichiometry()
    elif choice == "3":
      solution_concentration()
    elif choice == "4":
      unit_conversions()
    elif choice == "5":
      percent_yield()
    elif choice == "6":
      percentage_per_mass()
    elif choice == "7":
      exit_calculator()
    else:
      print("\ninvalid choice, please select a valid option")

if __name__ == "__main__":
  # Import Streamlit at the top level
  import streamlit as st

  def main_streamlit():
    st.title("Chemistry Calculator")
    st.write("Welcome to the Chemistry Calculator! Select an operation from the sidebar to get started.")
    
    menu_options = [
      "Mole & Mass Calculations", 
      "Stoichiometry", 
      "Solution Concentration", 
      "Unit Conversions", 
      "Percent Yield Calculation", 
      "Percentage Per Mass Calculation"
    ]

    option = st.sidebar.selectbox("Select an operation", menu_options)

    if option == "Mole & Mass Calculations":
      st.header("Mole & Mass Calculations")

      sub_option = st.selectbox("Choose calculation type:", 
                              ["Find moles", "Find mass", "Find compound molar mass", "Element molar mass"])

      if sub_option == "Find moles":
        st.subheader("Finding number of moles")
        mass = st.number_input("Mass (g):", placeholder="1.0", value=None, min_value=0.0,step=0.1)
        element = st.text_input("Element name or symbol:", placeholder="H").capitalize()

        if st.button("Calculate"):
          if element in elements:
            element_obj = elements[element]
            moles = mass/element_obj.atomic_mass
            st.success(f"Number of moles: {moles:.4f} mol")
          else:
            st.warning("Element not found. Please provide its molar mass.")
            molar_mass = st.number_input("Molar mass (g/mol):", min_value=0.0001, value=1.0, step=0.1)
            if st.button("Calculate with custom molar mass"):
                moles = mass/molar_mass
                st.success(f"Number of moles: {moles:.4f} mol")

      elif sub_option == "Find mass":
        st.subheader("Finding mass")
        moles = st.number_input("Number of moles (mol):", placeholder="1.0", value=None, min_value=0.0, step=0.1)
        element = st.text_input("Element name or symbol:", placeholder="H").capitalize()

        if st.button("Calculate"):
          if element in elements:
            element_obj = elements[element]
            mass = moles*element_obj.atomic_mass
            st.success(f"Mass: {mass:.4f} g")
          else:
            st.warning("Element not found. Please provide its molar mass.")
            molar_mass = st.number_input("Molar mass (g/mol):", min_value=0.0001, value=1.0, step=0.1)
            if st.button("Calculate with custom molar mass"):
                mass = moles*molar_mass
                st.success(f"Mass: {mass:.4f} g")

      elif sub_option == "Find compound molar mass":
        st.subheader("Compound Molar Mass")

        if 'compound_mass' not in st.session_state:
            st.session_state.compound_mass = 0
            st.session_state.elements_added = []

        st.write("Current compound molar mass: ", st.session_state.compound_mass, "g/mol")

        if st.session_state.elements_added:
            st.write("Elements added:")
            for elem in st.session_state.elements_added:
                st.write(f"- {elem}")

        st.write("Add elements to your compound:")
        col1, col2 = st.columns(2)

        with col1:
          element = st.text_input("Element name or symbol:", placeholder="C").capitalize()
        with col2:
          atoms = st.number_input("Number of atoms:", placeholder="1", value=None, min_value=1, step=1)

        if st.button("Add Element"):
          if element in elements:
            element_obj = elements[element]
            mass_contribution = atoms * element_obj.atomic_mass
            st.session_state.compound_mass += mass_contribution
            st.session_state.elements_added.append(f"{element} × {atoms} = {mass_contribution:.4f} g/mol")
            st.success(f"Added {element} × {atoms}. Current molar mass: {st.session_state.compound_mass:.4f} g/mol")
            st.rerun()
          else:
            st.error(f"{element} not found in the periodic table")

        if st.button("Reset Compound"):
            st.session_state.compound_mass = 0
            st.session_state.elements_added = []
            st.success("Compound reset")
            st.rerun()

      elif sub_option == "Element molar mass":
        st.subheader("Element Molar Mass")
        element = st.text_input("Element name or symbol:", placeholder="Fe").capitalize()

        if st.button("Search"):
          if element in elements:
            element_obj = elements[element]
            st.success(f"{element} molar mass: {element_obj.atomic_mass} g/mol")
          else:
            st.error(f"{element} does not correspond to any recognized chemical element")

    elif option == "Stoichiometry":
      st.header("Stoichiometry")

      sub_option = st.selectbox("Choose calculation type:", 
                              ["Empirical & Molecular Formulas", "Formulas of Hydrates"])
      # Improved empirical formula calculation for Streamlit
# To be added to your main_streamlit() function under the "Empirical & Molecular Formulas" section

      if sub_option == "Empirical & Molecular Formulas":
        st.subheader("Empirical & Molecular Formulas")

  # Initialize session state variables
        if 'elements_added' not in st.session_state:
          st.session_state.elements_added = []
    
  # Display current elements added with properly formatted information
        if st.session_state.elements_added:
          st.write("Elements added:")
    
    # Create a neat table to display elements
          element_data = []
          for elem, mass, moles in st.session_state.elements_added:
            element_data.append({"Element": elem, "Mass (g)": f"{mass:.4f}", "Moles": f"{moles:.4f}"})
    
          st.table(element_data)

  # Input fields for adding elements
        col1, col2 = st.columns(2)
        with col1:
          element = st.text_input("Element name or symbol:", "C").capitalize()
        with col2:
          mass = st.number_input("Mass of element (g):", min_value=0.01, value=10.0, step=0.01)

  # Add element to the list
        if st.button("Add Element"):
          if element in elements: 
            element_obj = elements[element]
            moles = mass / element_obj.atomic_mass
      # Store element, mass, and moles
            st.session_state.elements_added.append((element, mass, moles))
            st.success(f"{element} ({element_obj.name}) added: {mass:.4f}g = {moles:.4f} mol")
            st.rerun()
          else:
            st.error(f"{element} not found in the periodic table")

  # Calculate empirical formula
        if st.button("Calculate Empirical Formula"):
          if len(st.session_state.elements_added) >= 1:
      # Extract moles for all elements
            moles_dict = {elem: moles for elem, _, moles in st.session_state.elements_added}
      
      # Find the minimum moles
            min_moles = min(moles_dict.values())
      
      # Calculate the ratio for each element
            ratio_dict = {elem: moles/min_moles for elem, moles in moles_dict.items()}
      
      # Find smallest multiplier to get whole numbers
            decimals = [val % 1 for val in ratio_dict.values() if val % 1 != 0]
            if decimals:
              denominators = [round(1/d) for d in decimals]
              lcm = denominators[0]
              for num in denominators[1:]:
                lcm = (lcm*num)//gcd(lcm, num)
              multiplier = lcm
            else:
              multiplier = 1
      
      # Calculate final whole-number ratios
            final_ratios = {elem: round(num*multiplier) for elem, num in ratio_dict.items()}
      
      # Generate the empirical formula with subscripts
            subscript_digits = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
            empirical_formula = "".join([f"{elem}{str(num).translate(subscript_digits) if num>1 else''}" 
                                  for elem, num in final_ratios.items()])
      
      # Display formula and ratio details
            st.success(f"Empirical formula: {empirical_formula}")
      
      # Display detailed ratio calculation
            st.write("### Calculation Details")
            details = []
            for elem, _, moles in st.session_state.elements_added:
              ratio = ratio_dict[elem]
              final = final_ratios[elem]
              details.append({
                "Element": elem,
                "Moles": f"{moles:.4f}",
                "Relative Ratio": f"{ratio:.2f}",
                "Final Ratio": final
        })
            st.table(details)
      
          else:
            st.error("Please add at least one element to calculate the empirical formula.")

        if st.button("Reset all elements"):
          st.session_state.elements_added = []
          st.success("All elements reset!")
          st.rerun()
        
      
      elif sub_option == "Formulas of Hydrates":
        st.write("coming soon")

      
    elif option == "Solution Concentration":
      st.header("Solution Concentration")

      sub_option = st.selectbox("Choose concentration type:", 
                              ["Molarity", "Molality", "Mole Fraction", "Mass Percentage", "Volume Percentage"])

      if sub_option == "Molarity":
        st.subheader("Molarity Calculation")
        calc_type = st.radio("Calculate:", ["Find molarity", "Find moles", "Find volume"])

        if calc_type == "Find molarity":
          moles = st.number_input("Solute moles (mol):", min_value=0.0, value=1.0, step=0.1)
          volume = st.number_input("Solution volume (L):", min_value=0.0001, value=1.0, step=0.1)

          if st.button("Calculate"):
            molarity = moles/volume
            st.success(f"Molarity: {molarity:.4f} M")

        elif calc_type == "Find moles":
          molarity = st.number_input("Molarity (M):", min_value=0.0, value=1.0, step=0.1)
          volume = st.number_input("Solution volume (L):", min_value=0.0001, value=1.0, step=0.1)

          if st.button("Calculate"):
            moles = molarity * volume
            st.success(f"Solute moles: {moles:.4f} mol")

        elif calc_type == "Find volume":
          molarity = st.number_input("Molarity (M):", min_value=0.0001, value=1.0, step=0.1)
          moles = st.number_input("Solute moles (mol):", min_value=0.0, value=1.0, step=0.1)

          if st.button("Calculate"):
            volume = moles / molarity
            st.success(f"Solution volume: {volume:.4f} L")

      elif sub_option == "Molality":
        st.subheader("Molality Calculation")
        calc_type = st.radio("Calculate:", ["Find molality", "Find moles", "Find mass"])

        if calc_type == "Find molality":
          moles = st.number_input("Solute moles (mol):", min_value=0.0, value=1.0, step=0.1)
          mass = st.number_input("Solvent mass (kg):", min_value=0.0001, value=1.0, step=0.1)

          if st.button("Calculate"):
            molality = moles/mass
            st.success(f"Molality: {molality:.4f} mol/kg")

        elif calc_type == "Find moles":
          molality = st.number_input("Molality (mol/kg):", min_value=0.0, value=1.0, step=0.1)
          mass = st.number_input("Solvent mass (kg):", min_value=0.0001, value=1.0, step=0.1)

          if st.button("Calculate"):
            moles = molality * mass
            st.success(f"Solute moles: {moles:.4f} mol")

        elif calc_type == "Find mass":
          molality = st.number_input("Molality (mol/kg):", min_value=0.0001, value=1.0, step=0.1)
          moles = st.number_input("Solute moles (mol):", min_value=0.0, value=1.0, step=0.1)

          if st.button("Calculate"):
            mass = moles / molality
            st.success(f"Solvent mass: {mass:.4f} kg")

      elif sub_option == "Mole Fraction":
        st.subheader("Mole Fraction Calculation")
        moles_solute = st.number_input("Solute moles (mol):", min_value=0.0, value=1.0, step=0.1)
        moles_total = st.number_input("Total moles in solution (mol):", min_value=0.0001, value=10.0, step=0.1)

        if st.button("Calculate"):
          mole_fraction = moles_solute / moles_total
          st.success(f"Mole fraction: {mole_fraction:.4f}")

      elif sub_option == "Mass Percentage":
        st.subheader("Mass Percentage Calculation")
        mass_solute = st.number_input("Solute mass (g):", min_value=0.0, value=1.0, step=0.1)
        mass_solution = st.number_input("Solution mass (g):", min_value=0.0001, value=10.0, step=0.1)

        if st.button("Calculate"):
          mass_percentage = (mass_solute / mass_solution) * 100
          st.success(f"Mass percentage: {mass_percentage:.2f}%")

      elif sub_option == "Volume Percentage":
        st.subheader("Volume Percentage Calculation")
        volume_solute = st.number_input("Solute volume (mL):", min_value=0.0, value=1.0, step=0.1)
        volume_solution = st.number_input("Solution volume (mL):", min_value=0.0001, value=10.0, step=0.1)

        if st.button("Calculate"):
          volume_percentage = (volume_solute / volume_solution) * 100
          st.success(f"Volume percentage: {volume_percentage:.2f}%")

    elif option == "Unit Conversions":
      st.header("Unit Conversions")

      sub_option = st.selectbox("Choose conversion type:", 
                               ["Moles-Particles", "Liters-Moles", "Temperature", "Energy", "Pressure"])

      if sub_option == "Moles-Particles":
        st.subheader("Moles-Particles Conversion")
        conv_type = st.radio("Convert:", ["Moles to Particles", "Particles to Moles"])

        if conv_type == "Moles to Particles":
          moles = st.number_input("Moles (mol):", min_value=0.0, value=1.0, step=0.1)

          if st.button("Convert"):
            avogadro = 6.0221367e23
            particles = moles * avogadro
            st.success(f"Number of particles: {particles:.4e}")

        elif conv_type == "Particles to Moles":
          particles = st.number_input("Number of particles:", min_value=0.0, value=6.02e23, step=1.0e22, format="%.2e")

          if st.button("Convert"):
            avogadro = 6.0221367e23
            moles = particles / avogadro
            st.success(f"Number of moles: {moles:.4f} mol")

      elif sub_option == "Temperature":
        st.subheader("Temperature Conversion")
        conv_type = st.radio("Convert:", ["Kelvin to Celsius", "Celsius to Kelvin"])

        if conv_type == "Kelvin to Celsius":
          kelvin = st.number_input("Temperature in Kelvin:", value=298.15, step=0.1)

          if st.button("Convert"):
            celsius = kelvin - 273.15
            st.success(f"Temperature in Celsius: {celsius:.2f} °C")

        elif conv_type == "Celsius to Kelvin":
          celsius = st.number_input("Temperature in Celsius:", value=25.0, step=0.1)

          if st.button("Convert"):
            kelvin = celsius + 273.15
            st.success(f"Temperature in Kelvin: {kelvin:.2f} K")

      elif sub_option == "Energy":
        st.subheader("Energy Conversion")
        conv_type = st.radio("Convert:", ["Joules to Calories", "Calories to Joules"])

        if conv_type == "Joules to Calories":
          joules = st.number_input("Energy in Joules (J):", value=4.184, step=0.1)

          if st.button("Convert"):
            calories = joules / 4.184
            st.success(f"Energy in Calories: {calories:.4f} Cal")

        elif conv_type == "Calories to Joules":
          calories = st.number_input("Energy in Calories (Cal):", value=1.0, step=0.1)

          if st.button("Convert"):
            joules = calories * 4.184
            st.success(f"Energy in Joules: {joules:.4f} J")

      elif sub_option == "Pressure":
        st.subheader("Pressure Conversion")
        conv_type = st.selectbox("Convert:", ["atm to Pa", "Pa to atm", "mmHg to Pa", "Pa to mmHg", "atm to mmHg", "mmHg to atm"])

        if "atm to Pa" in conv_type:
          pressure = st.number_input(f"Pressure in atm:", value=1.0, step=0.1)
          if st.button("Convert"):
            result = pressure * 101325
            st.success(f"Pressure in Pa: {result:.2f} Pa")
        elif "Pa to atm" in conv_type:
          pressure = st.number_input(f"Pressure in Pa:", value=101325.0, step=1000.0)
          if st.button("Convert"):
            result = pressure / 101325
            st.success(f"Pressure in atm: {result:.6f} atm")
        elif "mmHg to Pa" in conv_type:
          pressure = st.number_input(f"Pressure in mmHg:", value=760.0, step=10.0)
          if st.button("Convert"):
            result = pressure * 133.322
            st.success(f"Pressure in Pa: {result:.2f} Pa")
        elif "Pa to mmHg" in conv_type:
          pressure = st.number_input(f"Pressure in Pa:", value=101325.0, step=1000.0)
          if st.button("Convert"):
            result = pressure / 133.322
            st.success(f"Pressure in mmHg: {result:.2f} mmHg")
        elif "atm to mmHg" in conv_type:
          pressure = st.number_input(f"Pressure in atm:", value=1.0, step=0.1)
          if st.button("Convert"):
            result = pressure * 760
            st.success(f"Pressure in mmHg: {result:.2f} mmHg")
        elif "mmHg to atm" in conv_type:
          pressure = st.number_input(f"Pressure in mmHg:", value=760.0, step=10.0)
          if st.button("Convert"):
            result = pressure / 760
            st.success(f"Pressure in atm: {result:.6f} atm")

    elif option == "Percent Yield Calculation":
      st.header("Percent Yield Calculation")

      actual_yield = st.number_input("Actual yield (g):", min_value=0.0, value=1.0, step=0.1)
      theoretical_yield = st.number_input("Theoretical yield (g):", min_value=0.0001, value=1.0, step=0.1)

      if st.button("Calculate"):
        percentage = actual_yield/theoretical_yield*100
        st.success(f"Percent yield: {percentage:.2f}%")

    elif option == "Percentage Per Mass Calculation":
      st.header("Percentage Per Mass Calculation")

      element_mass = st.number_input("Element mass (g):", min_value=0.0, value=1.0, step=0.1)
      compound_mass = st.number_input("Compound mass (g):", min_value=0.0001, value=10.0, step=0.1)

      if st.button("Calculate"):
        percentage = element_mass*100/compound_mass
        st.success(f"Percentage per mass: {percentage:.2f}%")

    st.markdown("\n\n\n\n \nDeveloped by Omamah Fath AL-Rahman Alsunni")


  # Run the Streamlit app
  main_streamlit()
