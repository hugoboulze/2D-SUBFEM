#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create .mat file (format Zset/Zebulon) for three different rheologies:
    - elastic
    - maxwell
    - burgers
"""

import os

# ---------------------------------------------------------------------------
# 2. FONCTIONS UTILITAIRES
# ---------------------------------------------------------------------------

def format_table(table, width_val=10, width_depth=12):
    """Formate une table (valeur, profondeur) en lignes alignees."""
    lines = []
    for val, depth in table:
        lines.append(f"{val:<{width_val}.4g}  {depth:<{width_depth}.6g}")
    return "\n".join(lines)


def build_C_table(G_table, S):
    """C = G * 3 * S, meme profondeurs que G_table."""
    return [(g_val * 3.0 * S, depth) for g_val, depth in G_table]

MAT_CHOICES = ("elastic", "maxwell", "burgers")

# ---------------------------------------------------------------------------
# 3. CONSTRUCTION DU CONTENU DU FICHIER
# ---------------------------------------------------------------------------

def write_mat_file(mat, eta_M=None, G_table=None, K_table=None, eta_K=None, S=1, elastic_mode="layered",
                    young=None, poisson=None, output_path="asthenosphere.mat"):
    """
    Write .mat file (format Zset).

    Parameters
    ----------
    mat : str
        "elastic", "maxwell" ou "burgers".
    eta_M : float
        Viscosite de Maxwell [Pa.s], utilisee dans **potential gen_evp ep.
    eta_K : float
        Viscosite de Kelvin-Voigt [Pa.s], utilisee dans **potential gen_evp ei
        (utilisee uniquement si mat == "burgers").
    S : float
        Facteur scalaire tel que C = G * 3 * S (table cinematique du potentiel ei).
    elastic_mode : str
        "layered" (par defaut) : elasticite donnee par les tables G_table/K_table
        en fonction de la profondeur.
        "constant" : elasticite donnee par un seul module de Young et un seul
        coefficient de Poisson (necessite young et poisson).
    young : float, optionnel
        Module de Young [Pa], requis si elastic_mode == "constant".
    poisson : float, optionnel
        Coefficient de Poisson [-], requis si elastic_mode == "constant".
    output_path : str
        Chemin du fichier .mat a ecrire.
    """
    if mat not in MAT_CHOICES:
        raise ValueError(f"mat must be on of these: {MAT_CHOICES}, received : {mat!r}")
    if elastic_mode not in ("layered", "constant"):
        raise ValueError(f"elastic_mode doit etre 'layered' ou 'constant', recu : {elastic_mode!r}")
    if elastic_mode == "constant" and (young is None or poisson is None):
        raise ValueError("Young modulus and Poisson ratio are required in 'constant' mode")

    # -- en-tete + elasticite (toujours presents) --------------------------
    
    parts = [f"""%\n***behavior gen_evp"""]
    
    if elastic_mode == "layered":
        
        elasticity_block = f"""\n**elasticity isotropic
\n#shear modulus G [Pa] 
G       profondeur_m
{format_table(G_table)}
\n#bulk modulus K [Pa] 
K           profondeur_m
{format_table(K_table)}"""
        
        G_scalar = None
    
    else:
        elasticity_block = f"""\n**elasticity isotropic
  young {young:.6g}
  poisson {poisson:.6g}"""
   
        G_scalar = young / (2.0 * (1.0 + poisson))

    parts.append(elasticity_block)
    if mat!='elastic':
        # -- potentiel gen_evp ep (Maxwell element) -------------------------------------
        if mat in ("maxwell", "burgers"):
            
            parts.append(f"""\n% ------------------- %
    %   MAXWELL (dashpot)
    % ------------------- %""")
            
            parts.append(f"""\n**potential gen_evp ep
     *flow norton
       n 1
       K  function (3*{eta_M:.6g})/86400 ;
     *isotropic constant
       R0 100""")
    
        # -- potentiel gen_evp ei (Kelvin-Voigt element) ------------------
        if mat == "burgers":
            
            if elastic_mode == "layered":
                C_block = f"""C    profondeur_m
    {format_table(build_C_table(G_table, S))}"""
            
            else:
                C_value = G_scalar * 3.0 * S
                C_block = f"C    {C_value:.6g}"
            
            parts.append(f"""\n% ------------------- %
    %   KELVIN-VOIGT body
    % ------------------- %""")
            
            parts.append(f"""\n**potential gen_evp ei
     *criterion mises
     *flow norton
       n  1
       K  function (3*{eta_K:.6g})/86400 ;
     *kinematic linear  xl
       {C_block}
     *isotropic constant
       R0 0.0""")

    parts.append("\n***return")

    content = "\n".join(parts)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)

    print(f"Fichier ecrit : {os.path.abspath(output_path)}")
    return content


# ---------------------------------------------------------------------------
# 4. EXEMPLE D'UTILISATION
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    
    # =============================================================================
    #  PREM tables
    # =============================================================================
    
    # --- Shear modulus table: G [Pa] vs profondeur [m] ---
    G_PREM = [
        (2.94e11, -3.00e6),  
        (2.94e11, -2.89e6),
        (1.93e11, -1.25e6),
        (1.55e11, -6.71e5),
        (1.25e11, -6.60e5),
        (9.06e10, -4.06e5),
        (8.24e10, -4.01e5),
        (6.56e10, -2.21e5),
        (68.e9,   -60e3),
        (68.e9,   -30.1e3),
        (26.e9,   -30e3),
        (26.e9,    1.e4),
    ]

    # --- Bulk modulus table: K [Pa] vs profondeur [m] ---
    K_PREM = [
        (6.56e11, -3.00e6),  
        (6.56e11, -2.89e6),
        (3.73e11, -1.25e6),
        (3.00e11, -6.71e5),
        (2.42e11, -6.60e5),
        (1.76e11, -4.06e5),
        (1.60e11, -4.01e5),
        (1.27e11, -2.21e5),
        (130.e9,  -60e3),
        (130.e9,  -30.1e3),
        (52.e9,   -30e3),
        (52.e9,    1.e4),
    ]
    
    # =============================================================================
    # EXAMPLE 1 - Elastic constant
    # =============================================================================
    write_mat_file(
        mat="elastic",
        elastic_mode='constant',
        young=50,
        poisson=0.25,
        output_path="./mat/elastic_constant.mat",
    )
    
    # =============================================================================
    # EXAMPLE 2 - Elastic layered
    # =============================================================================
    write_mat_file(
        mat="elastic",
        elastic_mode='layered',
        G_table=G_PREM,
        K_table=K_PREM,
        output_path="./mat/elastic_layered.mat",
    )
    
    # =============================================================================
    #  EXAMPLE 3 - Maxwell rheology
    # =============================================================================
    write_mat_file(
        mat="maxwell",
        elastic_mode="layered",
        G_table=G_PREM,
        K_table=K_PREM,
        eta_M=3.0e18,
        output_path="./mat/maxwell.mat",
    )
    
    # =============================================================================
    #  EXAMPLE 4 - Burgers rheology   
    # =============================================================================
    write_mat_file(
        mat="burgers",
        elastic_mode="layered",
        G_table=G_PREM,
        K_table=K_PREM,
        eta_M=3.0e19,
        eta_K=2.7e18,
        S=0.2, # mu_K (C) = S * mu_M (G)
        output_path="burgers.mat",
    )