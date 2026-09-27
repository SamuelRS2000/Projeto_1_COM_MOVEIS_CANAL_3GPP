import os
import sys

# Garante a importacao da subpasta 'Classes' no diretorio real deste script
_diretorio_raiz = os.path.dirname(os.path.abspath(__file__))
_pasta_classes = os.path.join(_diretorio_raiz, "Classes")

if _pasta_classes not in sys.path:
    sys.path.insert(0, _pasta_classes)

# Importa a ponta da cadeia de heranca
from calc_desv_doopler_fase import calc_desv_doopler
import matplotlib.pyplot as plt
import numpy as np


class Canal3GPP(calc_desv_doopler):
    """Classe consolidadora do Canal 3GPP TR 38.901.

    Concentra por heranca todas as caracteristicas:
      - Geometria e posicionamento BS/UE (gerar_param)
      - Parametros de Larga Escala e LoS/NLoS (gerar_param_larga_escala)
      - Atrasos multipercurso e potencias (gerar_atrasos)
      - Direcoes de chegada e vetores r_n (gerar_direcoes_chegadas)
      - Desvio Doppler e fases variante no tempo (calc_desv_doopler)
    """

    def __init__(self, ambiente="UMi", freq=3.5, N_multi=100, **kwargs):
        super().__init__(
            ambiente=ambiente, freq=freq, N_multi=N_multi, **kwargs
        )

    # Metodo utilitario para imprimir o resumo completo das propriedades do canal
    def exibir_resumo(self):
        print(f"\n{'='*55}")
        print(
            f"          RESUMO DO CANAL 3GPP TR 38.901 ({self.ambiente})          "
        )
        print(f"{'='*55}")
        print(
            f"Condicao de Enlace:       {'LoS (Visada Direta)' if self.flag_Loss == 1 else 'NLoS'}"
        )
        print(
            f"Frequencia Central:       {self.freq:.2f} GHz (lambda = {3e8/(self.freq*1e9):.4f} m)"
        )
        print(
            f"Posicao Relativa UE:      x={self.x:.1f} m, y={self.y:.1f} m, z={self.z:.1f} m"
        )
        print(
            f"Distancia 2D / 3D:        {self.d_BS_UE:.2f} m / {self.d_3D:.2f} m"
        )
        print(
            f"Velocidade do Terminal:   {self.mod_v_tx:.1f} km/h ({self.mod_v_tx/3.6:.2f} m/s)"
        )
        print(f"Vetor Direcao Velocidade: {self.vet_v_tx}")
        print(f"{'-'*55}")
        print(f"Espalhamento Atraso (DS): {self.espa_atraso * 1e9:.2f} ns")
        print(f"Fator de Rice (K):        {self.rice_factor:.2f}")
        print(f"Numero de Multipercursos: {self.N_multi}")
        print(
            f"Intervalo de Atrasos:     [{self.tau_n[0]*1e9:.2f}, {self.tau_n[-1]*1e9:.2f}] ns"
        )
        print(
            f"Desvio Doppler Min / Max: [{np.min(self.nu_n):.2f}, {np.max(self.nu_n):.2f}] Hz"
        )
        print(f"Soma Total das Potencias: {np.sum(self.alpha_sq):.4f}")
        print(f"{'='*55}\n")


# --- Teste de instanciacao unica do canal consolidado ---
if __name__ == "__main__":

    # Cria o objeto consolidado com todas as caracteristicas
    canal = Canal3GPP(
        ambiente="UMi",
        freq=3.5,
        N_multi=100,
        x_UE=30.0,
        y_UE=40.0,
        mod_v_tx=3.0,
    )

    # Exibe o sumario de todos os parametros gerados
    canal.exibir_resumo()

    # Todos os metodos de plotagem das classes base estao prontos para uso:
    canal.plot_PDP()  # Perfil de Atraso de Potencia (Slide 15)
    canal.plot_doppler()  # Espectro Doppler (Slide 26)
    canal.plot_r_n()  # Vetores de Incidencia 3D (Slide 24)