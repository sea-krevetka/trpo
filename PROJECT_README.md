# BLAS Custom Implementation with CI/CD

Этот проект демонстрирует:
- Реализацию пользовательской функции BLAS (TRMM - Triangular Matrix Multiplication)
- Интеграцию с OpenBLAS библиотекой
- Интерфейсные тесты с использованием CBLAS API
- Бенчмарки производительности
- Автоматическую интеграцию и развертывание (CI/CD) с GitHub Actions

## Структура проекта

```
trpo/
├── custom_func/
│   ├── custom_trmm.h      - Пользовательская реализация TRMM
│   ├── benchmark.cpp      - Бенчмарк для оценки производительности
│   └── Makefile
├── test_cblas_l3/
│   ├── test_cblas_level3.c - Интерфейсные тесты CBLAS L3
│   ├── fake_cblas.c        - Имитация BLAS функций для тестирования
│   ├── cblas.h             - CBLAS заголовок
│   ├── Makefile            - Правила сборки
│   └── README.md
├── OpenBLAS/               - Исходный код OpenBLAS библиотеки
└── .github/workflows/
    └── ci.yml              - GitHub Actions CI/CD конфигурация
```

## GitHub Actions CI/CD Pipeline

### Автоматические проверки

Pipeline автоматически запускается при каждом push в ветку `main` и выполняет:

1. **Установка зависимостей** - Необходимые библиотеки и компиляторы
2. **Сборка проекта** - Компиляция OpenBLAS и тестов
3. **Интерфейсные тесты** - Проверка корректности реализации
4. **Проверка пользовательской реализации** - Валидация custom_trmm.h
5. **Отчет о результатах** - Обновление статуса сборки

### Результаты CI/CD

- ✅ **Успешные сборки** - Все тесты пройдены
- ❌ **Неудачные сборки** - Обнаружены ошибки в тестах

## Как использовать

### Локальная сборка и тестирование

```bash
# Построить и запустить тесты
cd test_cblas_l3
make run_all

# Запустить бенчмарк
cd ../custom_func
g++ -O3 -march=native -fopenmp -o benchmark benchmark.cpp -lstdc++ -lm
./benchmark
```

### Размещение на GitHub

1. Создайте репозиторий на GitHub
2. Добавьте remote:
   ```bash
   git remote add origin https://github.com/YOUR_USERNAME/trpo.git
   git branch -M main
   git push -u origin main
   ```

3. GitHub Actions автоматически начнет выполнение при push

## Демонстрация CI/CD

### Успешный коммит (passing tests)
Обычный коммит с правильной реализацией - все тесты пройдены ✅

### Неудачный коммит (failing tests)
Коммит с изменением реализации, вызывающим ошибку в тестах ❌

## Требования

- GCC/G++ компилятор
- Make
- OpenBLAS библиотека
- CBLAS интерфейс

## Авторы

Демонстрационный проект для демонстрации навыков CI/CD с GitHub Actions.
