import os
import sys

# Garante a importacao a partir do diretorio real deste arquivo
_diretorio_real = os.path.dirname(os.path.abspath(__file__))
if _diretorio_real not in sys.path:
    sys.path.insert(0, _diretorio_real)

import matplotlib.pyplot as plt
from gerar_param import gerar_param
import numpy as np


class gerar_param_larga_escala(gerar_param):

    def __init__(self, ambiente="UMi", freq=3, **kwargs):
        # Heranca da classe base gerar_param
        super().__init__(ambiente=ambiente, freq=freq, **kwargs)

        if self.flag_Loss == 1:
            # Espalhamento de atraso (DS) - Tabela 7.5-6
            mu_lgDS = -0.24 * np.log10(1 + self.freq) - 7.14
            sigma_lgDS = 0.38
            lgDS = np.random.normal(mu_lgDS, sigma_lgDS)
            self.espa_atraso = 10**lgDS  # [s]

            # Fator de Rice (K) - Tabela 7.5-6
            mu_K = 9.0  # [dB]
            sigma_K = 5.0  # [dB]
            K_dB = np.random.normal(mu_K, sigma_K)
            self.rice_factor = 10 ** (K_dB / 10.0)  # Linear

            # Espalhamento Azimutal de Saída (ASD) - Tabela 7.5-6
            mu_lgASD = -0.05 * np.log10(1 + self.freq) + 1.21
            sigma_lgASD = 0.41
            lgASD = np.random.normal(mu_lgASD, sigma_lgASD)
            self.esp_Azt_out = 10**lgASD  # [graus]

            # Espalhamento Azimutal de Chegada (ASA) - Tabela 7.5-6
            mu_lgASA = -0.08 * np.log10(1 + self.freq) + 1.73
            sigma_lgASA = 0.014 * np.log10(1 + self.freq) + 0.28
            lgASA = np.random.normal(mu_lgASA, sigma_lgASA)
            self.esp_Azt_in = 10**lgASA  # [graus]

            # Espalhamento em Elevação de Chegada (ZSA/ESA) - Tabela 7.5-6
            mu_lgZSA = -0.17 * np.log10(1 + self.freq) + 0.73
            sigma_lgZSA = -0.01 * np.log10(1 + self.freq) + 0.34
            lgZSA = np.random.normal(mu_lgZSA, sigma_lgZSA)
            self.esp_elv_in = 10**lgZSA  # [graus]

            # Espalhamento em Elevação de Saída (ZSD/ESD) - Tabela 7.5-7
            mu_lgZSD = max(
                -0.5,
                -2.1 * (self.d_BS_UE / 1000.0)
                - 0.01 * (self.h_UT - 1.5)
                + 0.75,
            )
            sigma_lgZSD = 0.40
            lgZSD = np.random.normal(mu_lgZSD, sigma_lgZSD)
            self.esp_elv_out = 10**lgZSD  # [graus]

            # Sombreamento (Shadow Fading) - Tabela 7.5-6
            self.sf = np.random.normal(0, 4.0)  # [dB]

        else:
            # Espalhamento de atraso (DS) - Tabela 7.5-6
            mu_lgDS = -0.24 * np.log10(1 + self.freq) - 6.83
            sigma_lgDS = 0.16 * np.log10(1 + self.freq) + 0.28
            lgDS = np.random.normal(mu_lgDS, sigma_lgDS)
            self.espa_atraso = 10**lgDS  # [s]

            # Fator de Rice é nulo em NLoS
            self.rice_factor = 0.0

            # Espalhamento Azimutal de Saída (ASD) - Tabela 7.5-6
            mu_lgASD = -0.23 * np.log10(1 + self.freq) + 1.53
            sigma_lgASD = 0.11 * np.log10(1 + self.freq) + 0.33
            lgASD = np.random.normal(mu_lgASD, sigma_lgASD)
            self.esp_Azt_out = 10**lgASD  # [graus]

            # Espalhamento Azimutal de Chegada (ASA) - Tabela 7.5-6
            mu_lgASA = -0.08 * np.log10(1 + self.freq) + 1.81
            sigma_lgASA = 0.05 * np.log10(1 + self.freq) + 0.30
            lgASA = np.random.normal(mu_lgASA, sigma_lgASA)
            self.esp_Azt_in = 10**lgASA  # [graus]

            # Espalhamento em Elevação de Chegada (ZSA/ESA) - Tabela 7.5-6
            mu_lgZSA = -0.19 * np.log10(1 + self.freq) + 0.96
            sigma_lgZSA = 0.38
            lgZSA = np.random.normal(mu_lgZSA, sigma_lgZSA)
            self.esp_elv_in = 10**lgZSA  # [graus]

            # Espalhamento em Elevação de Saída (ZSD/ESD) - Tabela 7.5-7
            mu_lgZSD = max(
                -0.5,
                -2.1 * (self.d_BS_UE / 1000.0)
                - 0.01 * (self.h_UT - 1.5)
                + 0.90,
            )
            sigma_lgZSD = 0.49
            lgZSD = np.random.normal(mu_lgZSD, sigma_lgZSD)
            self.esp_elv_out = 10**lgZSD  # [graus]

            # Sombreamento (Shadow Fading) - Tabela 7.5-6
            self.sf = np.random.normal(0, 7.82)  # [dB]

        # Limitação dos espalhamentos angulares
        self.esp_Azt_out = float(min(self.esp_Azt_out, 104.0))
        self.esp_Azt_in = float(min(self.esp_Azt_in, 104.0))
        self.esp_elv_out = float(min(self.esp_elv_out, 52.0))
        self.esp_elv_in = float(min(self.esp_elv_in, 52.0))


# --- Teste isolado do arquivo gerar_param_larga_escala.py ---
if __name__ == "__main__":
    canal = gerar_param_larga_escala(
        ambiente="UMi", freq=3.5, x_UE=25.0, y_UE=30.0
    )
    print("\n--- Teste Larga Escala (UMi) ---")
    print("Ambiente:", canal.ambiente)
    print(f"Condição: {'LoS' if canal.flag_Loss == 1 else 'NLoS'}")
    print(f"Distância 2D: {canal.d_BS_UE:.2f} m | Distância 3D: {canal.d_3D:.2f} m")
    print(
        f"Espalhamento de Atraso (DS): {canal.espa_atraso * 1e9:.2f} ns | Fator Rice (K): {canal.rice_factor:.2f}"
    )
    print(
        f"Espalhamento Azimutal (Saída/Chegada): {canal.esp_Azt_out:.2f}° / {canal.esp_Azt_in:.2f}°"
    )
    print(
        f"Espalhamento Elevação (Saída/Chegada): {canal.esp_elv_out:.2f}° / {canal.esp_elv_in:.2f}°"
    )
    print(f"Sombreamento (SF): {canal.sf:.2f} dB")