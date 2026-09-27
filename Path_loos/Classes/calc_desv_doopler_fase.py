import os
import sys

# Garante a importacao a partir do diretorio real deste arquivo
_diretorio_real = os.path.dirname(os.path.abspath(__file__))
if _diretorio_real not in sys.path:
    sys.path.insert(0, _diretorio_real)

from gerar_direcoes_chegadas import gerar_direcoes_chegadas
import matplotlib.pyplot as plt
import numpy as np


class calc_desv_doopler(gerar_direcoes_chegadas):

    def __init__(self, ambiente="UMi", freq=3, N_multi=100, **kwargs):
        # Heranca da classe base gerar_direcoes_chegadas
        super().__init__(
            ambiente=ambiente, freq=freq, N_multi=N_multi, **kwargs
        )

        # Constantes de propagacao
        self.c = 3e8  # Velocidade da luz [m/s]
        self.f_c = self.freq * 1e9  # Frequencia central convertida para [Hz]
        self.lambda_c = self.c / self.f_c  # Comprimento de onda [m]

        # Calculo inicial do Doppler
        self.atualizar_velocidade(self.mod_v_tx)

    # Recalcula as propriedades variantes no tempo mantendo o canal espacial fixo
    def atualizar_velocidade(self, nova_velocidade_kmh):
        self.mod_v_tx = nova_velocidade_kmh
        
        # Velocidade escalar convertida de km/h para m/s
        v_rx_ms = self.mod_v_tx / 3.6

        # Desvio Doppler nu_n (Slide 25 - Eq. 9)
        # Produto interno entre cada vetor coluna de r_n (3 x N) e o vetor unitario de velocidade (3,)
        produto_escalar = np.dot(self.r_n.T, self.vet_v_tx)
        self.nu_n = (v_rx_ms / self.lambda_c) * produto_escalar  # [Hz]

        # Fases Multipercurso - Componente estatica (Slide 27 - Eq. 11)
        self.phi_bar_n = 2.0 * np.pi * (self.f_c + self.nu_n) * self.tau_n  # [rad]

    # Avaliacao da fase variante no tempo (Slide 27 - Eq. 10)
    def calcular_fases(self, t):
        t = np.atleast_1d(t)
        # Retorna matriz de dimensoes (N_multi, len(t))
        return self.phi_bar_n[:, None] - 2.0 * np.pi * self.nu_n[:, None] * t

    # Grafico do Espectro Doppler (Slide 26 / Figura 5)
    def plot_doppler(self, titulo=None, ax=None, limite_x=None):
        show_plot = False
        
        # Criacao da figura apenas se nenhum eixo externo for fornecido
        if ax is None:
            fig, ax = plt.subplots(figsize=(9, 5))
            show_plot = True

        # Plotagem de cada componente Doppler com a sua respectiva potencia multipercurso
        markerline, stemlines, baseline = ax.stem(
            self.nu_n, self.alpha_sq, linefmt="k-", markerfmt="k^", basefmt=" "
        )
        plt.setp(stemlines, "linewidth", 1.2)
        plt.setp(markerline, "markersize", 6)

        # Escala logaritmica no eixo vertical
        ax.set_yscale("log")
        y_min_plot = 1e-5
        ax.set_ylim([y_min_plot, max(1e-1, np.max(self.alpha_sq) * 1.5)])

        # Trava o eixo X globalmente se o parametro limite_x for passado
        if limite_x is not None:
            ax.set_xlim([-limite_x, limite_x])
        else:
            nu_lim = max(np.max(np.abs(self.nu_n)) * 1.15, 1.0)
            ax.set_xlim([-nu_lim, nu_lim])

        # Formatacao dos rotulos
        ax.set_xlabel(r"Domínio de Desvio Doppler - $\nu$ (Hz)", fontsize=12)
        ax.set_ylabel("Espectro Doppler", fontsize=12)
        
        # Definicao do titulo dinamicamente
        titulo_padrao = titulo if titulo else f"Espectro Doppler (v = {self.mod_v_tx} km/h)"
        ax.set_title(titulo_padrao, fontsize=12, fontweight="bold")
        
        ax.grid(True, which="both", linestyle="--", alpha=0.6)

        # Informacao do desvio Doppler maximo no quadro superior
        max_nu = np.max(np.abs(self.nu_n))
        ax.text(
            0.95,
            0.93,
            rf"$\nu_{{max}} \approx {max_nu:.1f}\text{{ Hz}}$",
            transform=ax.transAxes,
            fontsize=11,
            verticalalignment="top",
            horizontalalignment="right",
            bbox=dict(
                boxstyle="square,pad=0.3", facecolor="white", edgecolor="black"
            ),
        )

        # Exibicao do grafico apenas se foi criado internamente
        if show_plot:
            plt.tight_layout()
            plt.show()
            
        return ax


# --- Teste isolado do arquivo calc_desv_doopler.py ---
if __name__ == "__main__":
    
    # Definicao de tres cenarios de velocidade (Pedestre, Veiculo Urbano, Rodovia)
    velocidades_kmh = [3.0, 50.0, 120.0]

    # Criacao do painel unificado (3 linhas, 1 coluna)
    fig, eixos = plt.subplots(3, 1, figsize=(10, 12))

    print("\n--- Teste Desvio Doppler (UMi) - Painel Multivelocidade ---")
    
    # Instanciacao de um unico canal fisico para garantir que as potencias 
    # e angulos de chegada sejam exatamente os mesmos nos tres graficos.
    canal_base = calc_desv_doopler(ambiente="UMi", freq=3.5, N_multi=100)

    # 1. Determina a escala horizontal global com base na velocidade maxima
    canal_base.atualizar_velocidade(max(velocidades_kmh))
    limite_x_global = np.max(np.abs(canal_base.nu_n)) * 1.10

    # 2. Gera os graficos travando a escala X em limite_x_global
    for v, ax_atual in zip(velocidades_kmh, eixos):
        
        # Atualiza a velocidade do móvel mantendo o canal espacial estático
        canal_base.atualizar_velocidade(v)

        max_doppler = np.max(np.abs(canal_base.nu_n))
        print(
            f"Gerando Espectro - Velocidade: {v:5.1f} km/h | "
            f"Condição: {'LoS' if canal_base.flag_Loss == 1 else 'NLoS'} | "
            f"Max Doppler: {max_doppler:6.2f} Hz"
        )

        canal_base.plot_doppler(
            titulo=f"Espectro Doppler - Mobilidade: {v} km/h", 
            ax=ax_atual,
            limite_x=limite_x_global
        )

    # Ajuste manual das margens absolutas e do espacamento vertical (hspace)
    fig.subplots_adjust(top=0.95, bottom=0.08, left=0.1, right=0.95, hspace=0.45)
    
    plt.show()