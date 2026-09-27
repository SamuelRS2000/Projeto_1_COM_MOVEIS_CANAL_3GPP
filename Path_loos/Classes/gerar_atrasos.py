import os
import sys

# Garante a importacao a partir do diretorio real deste arquivo
_diretorio_real = os.path.dirname(os.path.abspath(__file__))
if _diretorio_real not in sys.path:
    sys.path.insert(0, _diretorio_real)

from gerar_param_larga_escala import gerar_param_larga_escala
import matplotlib.pyplot as plt
import numpy as np


class gerar_atrasos(gerar_param_larga_escala):

    def __init__(
        self,
        ambiente="UMi",
        freq=3,
        N_multi=100,
        override_espa_atraso=None,
        **kwargs
    ):
        # Heranca da classe base gerar_param_larga_escala
        super().__init__(ambiente=ambiente, freq=freq, **kwargs)

        self.N_multi = N_multi  # Número de componentes multipercurso

        # Sobrescrita manual do Espalhamento de Atraso (Delay Spread), se fornecido - (isso serve para plots o padrão e ser calculado)
        if override_espa_atraso is not None:
            self.espa_atraso = float(override_espa_atraso)

        # Fator de proporcionalidade r_tau (Tabela 7.5-6)
        if self.flag_Loss == 1:
            self.r_tau = 2.7
        else:
            self.r_tau = 3.0

        # Média da distribuição dos atrasos (Slide 12)
        self.mu_tau = self.r_tau * self.espa_atraso

        # Geração dos atrasos exponenciais e ordenação (Slide 12)
        tau_raw = np.random.exponential(scale=self.mu_tau, size=self.N_multi)
        tau_norm = tau_raw - np.min(tau_raw)
        self.tau_n = np.sort(tau_norm)

        # Potência Multipercurso (Slide 13 e 14)
        # Termo de sombreamento por cluster xi_n (Tabela 7.5-6: sigma_xi = 3 dB para UMi)
        self.sigma_xi = 3.0  # [dB]
        self.ksi = np.random.normal(0.0, self.sigma_xi, size=self.N_multi)

        # Potência preliminar de cada componente (Slide 13 - Eq. 2)
        termo_exp = -self.tau_n * (
            (self.r_tau - 1.0) / (self.r_tau * self.espa_atraso)
        )
        alpha_hat_sq = np.exp(termo_exp) * (10.0 ** (-self.ksi / 10.0))

        # Normalização das potências (Slide 14)
        self.alpha_sq = np.zeros(self.N_multi)

        if self.flag_Loss == 1:
            # Caso LoS: componente principal separada e componentes dispersas normalizadas
            omega_c = np.sum(alpha_hat_sq[1:])
            self.alpha_sq[1:] = (1.0 / (self.rice_factor + 1.0)) * (
                alpha_hat_sq[1:] / omega_c
            )
            self.alpha_sq[0] = self.rice_factor / (self.rice_factor + 1.0)
        else:
            # Caso NLoS: todas as componentes normalizadas pela soma total
            omega_c = np.sum(alpha_hat_sq)
            self.alpha_sq = alpha_hat_sq / omega_c

    # Gráfico do Perfil de Atraso de Potência (PDP)
    def plot_PDP(self, titulo=None, ax=None):
        show_plot = False
        
        # Cria a figura apenas se nenhum eixo externo for fornecido
        if ax is None:
            fig, ax = plt.subplots(figsize=(9, 5))
            show_plot = True

        tau_us = self.tau_n * 1e6  # Converte atrasos para microssegundos (us)

        markerline, stemlines, baseline = ax.stem(
            tau_us, self.alpha_sq, linefmt="k-", markerfmt="k^", basefmt=" "
        )
        plt.setp(stemlines, "linewidth", 1.2)
        plt.setp(markerline, "markersize", 6)

        # Definição do limite inferior de visualização do eixo Y
        y_min_plot = 1e-5
        
        ax.set_yscale("log")
        ax.set_ylim([y_min_plot, max(1e-1, np.max(self.alpha_sq) * 1.5)])
        
        # Filtra apenas os atrasos cujas potências são visíveis no gráfico
        taus_visiveis = tau_us[self.alpha_sq >= y_min_plot]
        if len(taus_visiveis) > 0:
            max_tau_visivel = np.max(taus_visiveis)
        else:
            max_tau_visivel = np.max(tau_us)
            
        # Define o eixo X com margem de -2% e estendendo 20% após o último atraso visível
        ax.set_xlim([-0.02 * max_tau_visivel, max_tau_visivel * 1.20])

        ax.set_xlabel(r"Domínio de Atraso - $\tau$ ($\mu$s)", fontsize=12)
        ax.set_ylabel("PDP", fontsize=12)
        
        # Define o título dinamicamente
        if titulo:
            ax.set_title(titulo, fontsize=12, fontweight="bold")
            
        ax.grid(True, which="both", linestyle="--", alpha=0.6)

        sigma_tau_ns = self.espa_atraso * 1e9
        ax.text(
            0.95,
            0.93,
            rf"$\sigma_\tau = {sigma_tau_ns:.2f}\text{{ ns}}$",
            transform=ax.transAxes,
            fontsize=12,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(
                boxstyle="square,pad=0.4", facecolor="white", edgecolor="black"
            ),
        )

        # Exibe o gráfico apenas se foi criado internamente
        if show_plot:
            plt.tight_layout()
            plt.show()
            
        return ax


# --- Teste isolado do arquivo gerar_atrasos.py ---
if __name__ == "__main__":
    
    # Definição dos três valores de espalhamento em nanosegundos
    valores_ds_ns = [10, 100, 1000]

    # Criação do painel unificado
    fig = plt.figure(figsize=(12, 8))
    
    # Estruturação dos eixos (2 colunas x 2 linhas)
    ax1 = plt.subplot(2, 2, 1)        # Canto superior esquerdo
    ax2 = plt.subplot(2, 2, 2)        # Canto superior direito
    ax3 = plt.subplot(2, 2, (3, 4))   # Linha inferior completa (ocupa colunas 1 e 2)
    
    eixos = [ax1, ax2, ax3]

    print("\n--- Teste Pequena Escala (UMi) - Painel Multigráfico ---")

    for ds_ns, ax_atual in zip(valores_ds_ns, eixos):
        ds_segundos = ds_ns * 1e-9
        
        canal_pequena = gerar_atrasos(
            ambiente="UMi", 
            freq=3.5, 
            N_multi=100, 
            override_espa_atraso=ds_segundos
        )

        print(f"Gerando PDP - Espalhamento: {ds_ns} ns | Condição: {'LoS' if canal_pequena.flag_Loss == 1 else 'NLoS'}")

        # Plota diretamente no eixo correspondente
        canal_pequena.plot_PDP(titulo=f"PDP com Espalhamento de Atraso: {ds_ns} ns", ax=ax_atual)

    # Ajusta o espaçamento entre os subplots e exibe o painel
    plt.tight_layout()
    plt.show()