#include <iostream>
#include <thread>
using std::thread;

void func1() {
    for (int i = 0; i < 3; i++) {
        std::cout << "쓰레드 1 작동중! \n";
    }
}

void func2() {
    for (int i = 0; i < 3; i++) {
        std::cout << "쓰레드 2 작동중! \n";
    }
}

void func3() {
    for (int i = 0; i < 3; i++) {
        std::cout << "쓰레드 3 작동중! \n";
    }
}
int main() {
    thread t1(func1);
    thread t2(func2);
    thread t3(func3);

    // t1.detach();
    // t2.detach();
    // t3.detach();

    t1.join();
    t2.join();
    t3.join();

    std::cout << "메인 함수 종료 \n";
}