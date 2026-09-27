import os
import sys

# Caminhos dinamicos relativos
_diretorio_raiz = os.path.dirname(os.path.abspath(__file__))
_pasta_path_loss = os.path.join(_diretorio_raiz, "Path_loos")
_pasta_classes = os.path.join(_diretorio_raiz, "Classes")

for pasta in [_diretorio_raiz, _pasta_path_loss, _pasta_classes]:
    if os.path.exists(pasta) and pasta not in sys.path:
        sys.path.insert(0, pasta)

from path_loss_3GPP_TR_38_901 import Canal3GPP
import matplotlib.pyplot as plt
import numpy as np


class CoerenciaCanal(Canal3GPP):

    def __init__(self, ambiente="UMi", freq=3.5, N_multi=100, override_espa_atraso=None, **kwargs):
        super().__init__(
            ambiente=ambiente, freq=freq, N_multi=N_multi, 
            override_espa_atraso=override_espa_atraso, **kwargs
        )

    # Funcao de Autocorrelacao Bidimensional Normalizada (Slide 37 - Eq. 13)
    def autocorrelacao_2D(self, kappa, sigma):
        kappa = np.atleast_1d(kappa)
        sigma = np.atleast_1d(sigma)
        omega_c = np.sum(self.alpha_sq)

        # Matrizes 3D para o Broadcasting (N, K, S)
        termo_freq = -2j * np.pi * self.tau_n[:, None, None] * kappa[None, :, None]
        termo_tempo = 2j * np.pi * self.nu_n[:, None, None] * sigma[None, None, :]

        soma = np.sum(
            self.alpha_sq[:, None, None] * np.exp(termo_freq + termo_tempo),
            axis=0,
        )
        return np.squeeze(soma / omega_c)

    # Recalculo limpo das potencias para variar o fator de Rice isoladamente
    def forcar_rice_factor(self, k_r_linear):
        self.flag_Loss = 1 if k_r_linear > 0 else 0
        self.rice_factor = float(k_r_linear)
        
        termo_exp = -self.tau_n * ((self.r_tau - 1.0) / (self.r_tau * self.espa_atraso))
        alpha_hat_sq = np.exp(termo_exp) * (10.0 ** (-self.ksi / 10.0))
        self.alpha_sq = np.zeros(self.N_multi)
        
        if self.flag_Loss == 1:
            omega_c = np.sum(alpha_hat_sq[1:])
            self.alpha_sq[1:] = (1.0 / (self.rice_factor + 1.0)) * (alpha_hat_sq[1:] / omega_c)
            self.alpha_sq[0] = self.rice_factor / (self.rice_factor + 1.0)
        else:
            omega_c = np.sum(alpha_hat_sq)
            self.alpha_sq = alpha_hat_sq / omega_c

    # Superficie 3D da Autocorrelacao
    def plot_autocorrelacao_3D(self, kappa_max=10e6, sigma_max=0.1, n_pontos=100):
        kappa = np.linspace(-kappa_max, kappa_max, n_pontos)
        sigma = np.linspace(-sigma_max, sigma_max, n_pontos)
        Z = np.abs(self.autocorrelacao_2D(kappa, sigma))
        K, S = np.meshgrid(kappa, sigma, indexing="ij")

        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(projection="3d")

        K_MHz = K * 1e-6
        S_ms = S * 1e3

        surf = ax.plot_surface(
            K_MHz, S_ms, Z, cmap="viridis", edgecolor="none", alpha=0.8, antialiased=True
        )

        ax.set_title("Função de Autocorrelação do Canal Bidimensional", fontsize=14, fontweight="bold", y=0.95)
        ax.set_xlabel(r"Desvio de Frequência $\kappa$ (MHz)", labelpad=10, fontsize=11)
        ax.set_ylabel(r"Desvio de Tempo $\sigma$ (ms)", labelpad=10, fontsize=11)
        ax.set_zlabel(r"$|\rho_{TT}(\kappa, \sigma)|$", labelpad=10, fontsize=11)
        
        ax.view_init(elev=25, azim=-45)
        fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.1, label="Correlação Absoluta")

        plt.tight_layout()
        plt.show()

    # --- Banda de Coerencia ---
    def calcular_banda_coerencia(self, n_pontos=3000):
        kappa = np.logspace(0, 10, n_pontos)
        omega_c = np.sum(self.alpha_sq)
        termo_exp = -2j * np.pi * self.tau_n[:, None] * kappa[None, :]
        rho_k0 = np.sum(self.alpha_sq[:, None] * np.exp(termo_exp), axis=0) / omega_c
        abs_rho_k0 = np.abs(rho_k0)

        bc_dict = {}
        for rho_b in [0.95, 0.8]:
            cruzou = np.where(abs_rho_k0 < rho_b)[0]
            bc_dict[rho_b] = kappa[cruzou[0]] if len(cruzou) > 0 else None
        return kappa, abs_rho_k0, bc_dict

    def plot_banda_coerencia(self, ax=None, titulo=None, n_pontos=3000):
        show_plot = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(9, 5))
            show_plot = True

        kappa, abs_rho, bc = self.calcular_banda_coerencia(n_pontos=n_pontos)

        ax.plot(kappa, abs_rho, color="#1f77b4", linewidth=2.0)
        ax.axhline(0.95, color="dimgray", linestyle="-.", linewidth=1.2)
        ax.axhline(0.80, color="dimgray", linestyle="-.", linewidth=1.2)

        titulo_bc = []
        if bc[0.95] is not None:
            ax.axvline(bc[0.95], color="dimgray", linestyle="--", linewidth=1.5)
            titulo_bc.append(rf"$B_C(0.95) = {bc[0.95]*1e-6:.1f}\text{{ MHz}}$")

        if bc[0.8] is not None:
            ax.axvline(bc[0.8], color="dimgray", linestyle="--", linewidth=1.5)
            titulo_bc.append(rf"$B_C(0.8) = {bc[0.8]*1e-6:.1f}\text{{ MHz}}$")

        ax.set_xscale("log")
        ax.set_xlim([1e0, 1e10])
        ax.set_ylim([0.2, 1.02])

        ax.set_title(titulo if titulo else ", ".join(titulo_bc), fontsize=12, fontweight="bold")
        ax.set_xlabel(r"Desvio de Frequência — $\kappa$ (Hz)", fontsize=11)
        ax.set_ylabel(r"$|\rho_{TT}(\kappa, 0)|$", fontsize=11)
        ax.grid(True, which="both", linestyle="--", alpha=0.5)

        info_str = "\n".join(titulo_bc)
        if info_str:
            ax.text(0.95, 0.93, info_str, transform=ax.transAxes, fontsize=11, 
                    va="top", ha="right", bbox=dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor="black"))

        if show_plot:
            plt.tight_layout()
            plt.show()

    # --- Tempo de Coerencia ---
    def calcular_tempo_coerencia(self, n_pontos=3000):
        sigma = np.logspace(-6, 0, n_pontos)
        omega_c = np.sum(self.alpha_sq)
        termo_exp = 2j * np.pi * self.nu_n[:, None] * sigma[None, :]
        rho_0s = np.sum(self.alpha_sq[:, None] * np.exp(termo_exp), axis=0) / omega_c
        abs_rho_0s = np.abs(rho_0s)

        tc_dict = {}
        for rho_t in [0.95, 0.8]:
            cruzou = np.where(abs_rho_0s < rho_t)[0]
            tc_dict[rho_t] = sigma[cruzou[0]] if len(cruzou) > 0 else None
        return sigma, abs_rho_0s, tc_dict

    def plot_tempo_coerencia(self, ax=None, titulo=None, n_pontos=3000):
        show_plot = False
        if ax is None:
            fig, ax = plt.subplots(figsize=(9, 5))
            show_plot = True

        sigma, abs_rho, tc = self.calcular_tempo_coerencia(n_pontos=n_pontos)

        ax.plot(sigma, abs_rho, color="#d62728", linewidth=2.0)
        ax.axhline(0.95, color="dimgray", linestyle="-.", linewidth=1.2)
        ax.axhline(0.80, color="dimgray", linestyle="-.", linewidth=1.2)

        titulo_tc = []
        if tc[0.95] is not None:
            ax.axvline(tc[0.95], color="dimgray", linestyle="--", linewidth=1.5)
            titulo_tc.append(rf"$T_C(0.95) = {tc[0.95]*1e3:.2f}\text{{ ms}}$")

        if tc[0.8] is not None:
            ax.axvline(tc[0.8], color="dimgray", linestyle="--", linewidth=1.5)
            titulo_tc.append(rf"$T_C(0.8) = {tc[0.8]*1e3:.2f}\text{{ ms}}$")

        ax.set_xscale("log")
        ax.set_xlim([1e-6, 1e0])
        ax.set_ylim([0.3, 1.02])

        ax.set_title(titulo if titulo else ", ".join(titulo_tc), fontsize=12, fontweight="bold")
        ax.set_xlabel(r"Desvio de Tempo — $\sigma$ (s)", fontsize=11)
        ax.set_ylabel(r"$|\rho_{TT}(0, \sigma)|$", fontsize=11)
        ax.grid(True, which="both", linestyle="--", alpha=0.5)

        info_str = "\n".join(titulo_tc)
        if info_str:
            ax.text(0.95, 0.93, info_str, transform=ax.transAxes, fontsize=11, 
                    va="top", ha="right", bbox=dict(boxstyle="square,pad=0.3", facecolor="white", edgecolor="black"))

        if show_plot:
            plt.tight_layout()
            plt.show()


# --- Execucao ---
if __name__ == "__main__":

    # =========================================================================
    # 0. SUPERFICIE 3D DA AUTOCORRELACAO 
    # =========================================================================
    print("\n--- Gerando Superficie 3D da Autocorrelacao (Feche a janela para continuar) ---")
    canal_3d = CoerenciaCanal(ambiente="UMi", freq=3.5, mod_v_tx=3.0, N_multi=100)
    canal_3d.plot_autocorrelacao_3D(kappa_max=10e6, sigma_max=0.1)

    # =========================================================================
    # 1. PAINEL 1 (3x1): BANDA DE COERENCIA VS ESPALHAMENTO DE ATRASO
    # =========================================================================
    print("\n--- Gerando Painel: Autocorrelação na Frequência vs Espalhamento de Atraso ---")
    fig1, axes1 = plt.subplots(3, 1, figsize=(9, 12))
    fig1.suptitle("Banda de Coerência vs Espalhamento de Atraso ($\sigma_\\tau$)", fontsize=15, fontweight="bold", y=0.96)
    
    ds_valores = [10, 100, 1000] # nanosegundos
    for i, ds_ns in enumerate(ds_valores):
        ds_segundos = ds_ns * 1e-9
        canal_ds = CoerenciaCanal(ambiente="UMi", freq=3.5, override_espa_atraso=ds_segundos)
        canal_ds.forcar_rice_factor(0) # Força NLoS para analisar atraso puro
        canal_ds.plot_banda_coerencia(ax=axes1[i], titulo=rf"Espalhamento de Atraso: $\sigma_\tau = {ds_ns}$ ns (NLoS)")
    
    # hspace aumentado de 0.45 para 0.65 para separar titulo do eixo X superior
    fig1.subplots_adjust(top=0.90, bottom=0.08, left=0.1, right=0.95, hspace=0.65)
    plt.show()

    # =========================================================================
    # 2. PAINEL 2 (3x1): BANDA DE COERENCIA VS FATOR DE RICE
    # =========================================================================
    print("\n--- Gerando Painel: Autocorrelação na Frequência vs Fator de Rice ---")
    fig2, axes2 = plt.subplots(3, 1, figsize=(9, 12))
    fig2.suptitle("Banda de Coerência vs Fator de Rice ($K_R$)", fontsize=15, fontweight="bold", y=0.96)
    
    rice_valores = [0, 10.0, 48.0] # Valores Lineares
    for i, kr in enumerate(rice_valores):
        canal_rice = CoerenciaCanal(ambiente="UMi", freq=3.5, override_espa_atraso=100e-9)
        canal_rice.forcar_rice_factor(kr)
        rotulo_rice = "NLoS" if kr == 0 else rf"$K_R = {kr}$"
        canal_rice.plot_banda_coerencia(ax=axes2[i], titulo=rf"Fator de Rice: {rotulo_rice} ($\sigma_\tau = 100$ ns)")

    # hspace aumentado para 0.65
    fig2.subplots_adjust(top=0.90, bottom=0.08, left=0.1, right=0.95, hspace=0.65)
    plt.show()

    # =========================================================================
    # 3. PAINEL 3 (3x1): TEMPO DE COERENCIA VS VELOCIDADE
    # =========================================================================
    print("\n--- Gerando Painel: Autocorrelação no Tempo vs Velocidade ---")
    fig3, axes3 = plt.subplots(3, 1, figsize=(9, 12))
    fig3.suptitle("Tempo de Coerência vs Mobilidade ($v_{rx}$)", fontsize=15, fontweight="bold", y=0.96)
    
    v_valores = [3.0, 50.0, 120.0] # km/h
    for i, v_kmh in enumerate(v_valores):
        canal_vel = CoerenciaCanal(ambiente="UMi", freq=3.5, mod_v_tx=v_kmh)
        canal_vel.plot_tempo_coerencia(ax=axes3[i], titulo=rf"Velocidade do Receptor: $v = {v_kmh}$ km/h")

    # hspace aumentado para 0.65
    fig3.subplots_adjust(top=0.90, bottom=0.08, left=0.1, right=0.95, hspace=0.65)
    plt.show()