# LDPC

Одна из задач спецкурса по теории кодирования. Программа оценивает сверху расстояние LDPC кода. 

Программа получает на вход весовую матрицу LDPC кода и вычисляет верхнюю оценку расстояния кода. Для этого в нескольких потоках вычисляются перманенты подматриц, а затем берется их минимальная сумма.
