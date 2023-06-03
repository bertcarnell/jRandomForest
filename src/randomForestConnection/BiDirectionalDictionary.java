package randomForestConnection;

import java.util.Dictionary;

public class BiDirectionalDictionary<TKey, TValue> : IDictionary<TKey, TValue>//, ICollection<KeyValuePair<TKey, TValue>>//, IEnumerable<KeyValuePair<TKey, TValue>>
{
    private Dictionary<TKey, TValue> oForwardDictionary;
    private Dictionary<TValue, TKey> oBackwardDictionary;

    /// <summary>
    /// 
    /// </summary>
    public BiDirectionalDictionary()
    {
        oForwardDictionary = new Dictionary<TKey, TValue>();
        oBackwardDictionary = new Dictionary<TValue, TKey>();
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="capacity"></param>
    public BiDirectionalDictionary(int capacity)
    {
        oForwardDictionary = new Dictionary<TKey, TValue>(capacity);
        oBackwardDictionary = new Dictionary<TValue, TKey>(capacity);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <param name="val"></param>
    public void Add(TKey key, TValue val)
    {
        oForwardDictionary.Add(key, val);
        oBackwardDictionary.Add(val, key);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <returns></returns>
    public TValue this[TKey key]
    {
        get
        {
            return oForwardDictionary[key];
        }
        set
        {
            // get the old value from the key to use on the backward dictionary
            TValue oldValue = oForwardDictionary[key];
            // set the value on the forward one
            oForwardDictionary[key] = value;
            // remove the instance of the old value from the backward dictionary (can't change the key at that location, so must remove)
            oBackwardDictionary.Remove(oldValue);
            // add the new pair
            oBackwardDictionary.Add(value, key);
        }
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="keys"></param>
    /// <returns></returns>
    public TValue[] getValueWithKey(TKey[] keys)
    {
        if (keys == null) return null;
        TValue[] results = new TValue[keys.Length];
        for (int i = 0; i < keys.Length; i++)
            results[i] = oForwardDictionary[keys[i]];
        return results;
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="val"></param>
    /// <returns></returns>
    public TKey getKeyWithValue(TValue val)
    {
        return oBackwardDictionary[val];
    }

    /// <summary>
    /// 
    /// </summary>
    /// <param name="vals"></param>
    /// <returns></returns>
    public TKey[] getKeyWithValue(TValue[] vals)
    {
        if (vals == null) return null;
        TKey[] results = new TKey[vals.Length];
        for (int i = 0; i < vals.Length; i++)
            results[i] = oBackwardDictionary[vals[i]];
        return results;
    }

    // Elements required by IDictionary

    /// <summary>
    /// 
    /// </summary>
    /// <param name="kvp"></param>
    public void Add(KeyValuePair<TKey, TValue> kvp)
    {
        oForwardDictionary.Add(kvp.Key, kvp.Value);
        oBackwardDictionary.Add(kvp.Value, kvp.Key);
    }
    /// <summary>
    /// 
    /// </summary>
    public void Clear()
    {
        oForwardDictionary.Clear();
        oBackwardDictionary.Clear();
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="kvp"></param>
    /// <returns></returns>
    public bool Contains(KeyValuePair<TKey, TValue> kvp)
    {
        return oForwardDictionary.Contains(kvp) & oBackwardDictionary.Contains(new KeyValuePair<TValue, TKey>(kvp.Value, kvp.Key));
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="kvp_array"></param>
    /// <param name="arrayIndex"></param>
    public void CopyTo(KeyValuePair<TKey, TValue>[] kvp_array, int arrayIndex)
    {
        foreach (KeyValuePair<TKey, TValue> kvp in oForwardDictionary)
        {
            kvp_array[arrayIndex] = kvp;
            arrayIndex++;
        }
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="kvp_array"></param>
    /// <param name="arrayIndex"></param>
    public void BackwardsCopyTo(KeyValuePair<TValue, TKey>[] kvp_array, int arrayIndex)
    {
        foreach (KeyValuePair<TValue, TKey> kvp in oBackwardDictionary)
        {
            kvp_array[arrayIndex] = kvp;
            arrayIndex++;
        }
    }
    /// <summary>
    /// 
    /// </summary>
    public int Count { get 
    { 
        Debug.Assert(oForwardDictionary.Count == oBackwardDictionary.Count);
        return oForwardDictionary.Count; 
    } }
    /// <summary>
    /// 
    /// </summary>
    public bool IsReadOnly { get { return false; } }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="kvp"></param>
    public bool Remove(KeyValuePair<TKey, TValue> kvp)
    {
        // single ampersand so both are done
        return oForwardDictionary.Remove(kvp.Key) & oBackwardDictionary.Remove(kvp.Value);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <returns></returns>
    public bool ContainsKey(TKey key)
    {
        return oForwardDictionary.ContainsKey(key);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <returns></returns>
    public bool BackwardContainsKey(TValue key)
    {
        return oBackwardDictionary.ContainsKey(key);
    }
    /// <summary>
    /// 
    /// </summary>
    public ICollection<TKey> Keys { get { return oForwardDictionary.Keys; } }
    /// <summary>
    /// 
    /// </summary>
    public ICollection<TValue> BackwardKeys { get { return oBackwardDictionary.Keys; } }
    /// <summary>
    /// 
    /// </summary>
    public ICollection<TValue> Values { get { return oForwardDictionary.Values; } }
    /// <summary>
    /// 
    /// </summary>
    public ICollection<TKey> BackwardValues { get { return oBackwardDictionary.Values; } }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <returns></returns>
    public bool Remove(TKey key)
    {
        // find the matching value
        TValue temp = oForwardDictionary[key];
        return oForwardDictionary.Remove(key) & oBackwardDictionary.Remove(temp);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public bool TryGetValue(TKey key, out TValue value)
    {
        return oForwardDictionary.TryGetValue(key, out value);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <param name="key"></param>
    /// <param name="value"></param>
    /// <returns></returns>
    public bool BackwardTryGetValue(TValue key, out TKey value)
    {
        return oBackwardDictionary.TryGetValue(key, out value);
    }
    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    public Dictionary<TKey, TValue>.Enumerator GetEnumerator()
    {
        return oForwardDictionary.GetEnumerator();
    }

    IEnumerator<KeyValuePair<TKey, TValue>> System.Collections.Generic.IEnumerable<System.Collections.Generic.KeyValuePair<TKey, TValue>>.GetEnumerator()
    {
        return oForwardDictionary.GetEnumerator();
    }

    IEnumerator System.Collections.IEnumerable.GetEnumerator()
    {
        return oForwardDictionary.GetEnumerator();
    }
    /// <summary>
    /// 
    /// </summary>
    /// <returns></returns>
    public Dictionary<TValue, TKey>.Enumerator GetBackwardEnumerator()
    {
        return oBackwardDictionary.GetEnumerator();
    }
}