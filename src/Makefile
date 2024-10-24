NAME = s21_matrix_oop.a

SRC = S21Matrix.cpp friends.cpp
HEADER = s21_matrix_oop.h
TEST = test.cpp

OBJ_FOLDER = obj
OBJ = $(SRC:%.cpp=$(OBJ_FOLDER)/%.o)

CLR = gcc
FLAGS = -std=c++17 -Werror -Wextra -Wall -I.

GCOV_REPORT_NAME = gcov_report

GREEN = \033[32m
YELLOW = \033[0;33m
RED = \033[0;31m
GRAY = \033[2;37m
CURSIVE = \033[3m
NO_COLOR = \x1b[0m


all: $(NAME)
.PHONY: all

$(NAME): $(OBJ_FOLDER) $(OBJ)
	@echo "\n"
	@ar rc $(NAME) obj/*	
	@ranlib $(NAME)
	@printf "$(CURSIVE)Assembling library ⚙️ ... $(YELLOW)%s$(NO_COLOR) successfully assembled! 🥂%-17s" $(NAME)""
	@echo "\n"

$(OBJ_FOLDER):
	@mkdir $(OBJ_FOLDER)

$(OBJ_FOLDER)/%.o: %.cpp $(HEADER) Makefile
	@$(CLR) -c $(FLAGS) $< -o $@
	@printf "$(CURSIVE)Compiling source files 🚀... $(YELLOW)%-13s $(NO_COLOR)---> %17.33s ✅\r" $< $@

clean:
	@echo "$(RED)Cleaning s21_matrix_oop.a objects..."
	@rm -rf obj 
	@echo "$(GRAY)$(CURSIVE)---s21_matrix_oop.a objects are cleaned!\n"
.PHONY: clean

clean_gcov:
	@echo "$(RED)Cleaning gcov / lcov staff..."
	@rm -rf $(GCOV_REPORT_NAME) $(GCOV_REPORT_NAME)_total.info $(GCOV_REPORT_NAME)_filtered.info *.gcda *.gcno report
	@echo "$(GRAY)$(CURSIVE)---gcov / lcov staff deleted!\n"
.PHONY: clean_gcov

fclean: clean clean_gcov
	@echo "$(RED)Cleaning $(NAME) & executable..."
	@rm -rf *.a a.out .*.swp *.*.gch .DS_Store test_bin a.out.dSYM $(NAME)
	@echo "$(GRAY)$(CURSIVE)---s21_matrix_oop.a staff is completely cleaned!\n"
.PHONY: fclean

re: fclean all
.PHONY: re

test: all
	@g++ -std=c++17 test.cpp $(NAME) -lgtest -lgtest_main -pthread -o test_bin
.PHONY: test

leaks: test
	@leaks --atExit -- ./test_bin
.PHONY: leaks

main: all
	@g++ -std=c++17 main.cpp $(NAME) && leaks --atExit -- ./a.out
.PHONY: main

correct_format:
	@clang-format -i $(SRC) $(TEST) $(HEADER)
.PHONY: correct_format

check_format:
	@clang-format --dry-run --Werror $(SRC) $(TEST) $(HEADER)
.PHONY: check_format

report/index.html: $(NAME)
	@g++ -std=c++17 -lgtest --coverage $(TEST) $(SRC) -o $(GCOV_REPORT_NAME)
	@./$(GCOV_REPORT_NAME)
	@lcov --no-external --no-recursion -c --directory . -o $(GCOV_REPORT_NAME)_total.info
	@genhtml -o report $(GCOV_REPORT_NAME)_total.info

open_gcov_report: report/index.html test.cpp
	@open report/index.html
.PHONY: open_gcov_report
