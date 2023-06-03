package randomForest;

import java.util.ArrayList;

public class ScootArray<T> {

    ArrayList<T> list;
    int pointer;
    
    public ScootArray() {
        list = new ArrayList<T>();
        pointer = 0;
    }
    
    public ScootArray(T[] data) {
        list = new ArrayList<T>();
        pointer = 0;
        for (int i = 0; i < data.length; i++) {
            list.add(data[i]);
        }
    }
    
    @SuppressWarnings("unchecked")
    public ScootArray(ScootArray<T> template) {
        list = (ArrayList<T>) template.list.clone();
        pointer = template.pointer;
    }
    
    public T get(int index) {
        return list.get(index + pointer);
    }
    
    public void set(int index, T data) {
        list.set(index + pointer, data);
    }
    
    public int size() {
        return list.size();
    }
    
    public ScootArray<T> scootLeft() {
        pointer--;
        return this;
    }
    
    public ScootArray<T> scootLeft(int num) {
        pointer -= num;
        return this;
    }
    
    public ScootArray<T> scootRight() {
        pointer++;
        return this;
    }
    
    public ScootArray<T> scootRight(int num) {
        pointer += num;
        return this;
    }
    
    public ScootArray<T> copy() {
        return new ScootArray<T>(this);
    }
    
}
