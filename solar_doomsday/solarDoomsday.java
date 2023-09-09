import java.lang.Math.sqrt;

// Not running correctly (Python code submitted)
public class Solution {
    public static int[] solution(int area) {
        ArrayList<int> ans = new ArrayList<int>();
        while (area > 0) {
            int x = (int)(sqrt((double)(area)));
            x *= x;
            area -= x;
            ans.add(x);
        }
        return ans.toArray();
    }
}
