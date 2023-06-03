package randomForestConnection;

public class ArrayPointer<T>
{
    /// <summary>
    /// 
    /// </summary>
    public T[] m_ptr; // have to expose this as public to use with the unmanaged calls

    /// <summary>
    /// 
    /// </summary>
    public T Value
    {
        get { return m_ptr[0]; }
        set { m_ptr[0] = value; }
    }

    /// <summary>
    /// 
    /// </summary>
    public ArrayPointer()
    {
        m_ptr = new T[1];
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="x"></param>
    public ArrayPointer(T x)
    {
        m_ptr = new T[1];
        m_ptr[0] = x;
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="val"></param>
    /// <returns></returns>
    public static implicit operator ArrayPointer<T>(T val)
    {
        return new ArrayPointer<T>(val);
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="val"></param>
    /// <returns></returns>
    public static implicit operator T(ArrayPointer<T> val)
    {
        Debug.Assert(val != null);
        return val.m_ptr[0];
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="val"></param>
    /// <returns></returns>
    public static explicit operator double(ArrayPointer<T> val)
    {
        return Convert.ToDouble(val.m_ptr[0]);
    }
}